#include "parametrizer.hpp"
#include "config.hpp"
#include "dedge.hpp"
#include "field-math.hpp"
#include "flow.hpp"
#include "localsat.hpp"
#include "optimizer.hpp"
#include "subdivide.hpp"

#include "dset.hpp"

#include <Eigen/Sparse>
#include <fstream>
#include <list>
#include <map>
#include <queue>
#include <set>

void Parametrizer::ComputeIndexMap(int with_scale) {
    // build edge info
    auto& V = hierarchy.mV[0];
    auto& F = hierarchy.mF;
    auto& Q = hierarchy.mQ[0];
    auto& N = hierarchy.mN[0];
    auto& O = hierarchy.mO[0];

    // ComputeOrientationSingularities();

    BuildEdgeInfo();

    std::vector<int> sharpv(V.cols(), 0);
    std::vector<int> sharpE(F.cols() * 3, 0);
    DisajointTree tree(V.cols());
    for (int i = 0; i < edge_diff.size(); ++i) {
        if (edge_diff[i][0] == 0 && edge_diff[i][1] == 0)
            tree.Merge(edge_values[i].x, edge_values[i].y);
    }
    for (int i = 0; i < sharp_edges.size(); ++i) {
        if (sharp_edges[i] == 0)
            continue;
        int v1 = F(i%3, i/3);
        int v2 = F((i+1)%3, i/3);
        sharpv[tree.Parent(v1)] = 1;
        sharpv[tree.Parent(v2)] = 1;
    }
    for (int i = 0; i < face_edgeIds.size(); ++i) {
        for (int j = 0; j < 3; ++j) {
            if (edge_diff[face_edgeIds[i][j]].array().abs().sum() == 1) {
                if (sharpv[tree.Parent(F(j, i))] && sharpv[tree.Parent(F((j+1)%3,i))])
                    sharpE[i * 3 + j] = 1;
            }
        }
    }
    std::swap(sharpE, sharp_edges);
    allow_changes.resize(edge_diff.size() * 2, 1);
    for (int i = 0; i < sharp_edges.size(); ++i) {
        if (sharp_edges[i] == 0)
            continue;
        int e = face_edgeIds[i/3][i%3];
        for (int k = 0; k < 2; ++k) {
            if (edge_diff[e][k] == 0)
                allow_changes[e * 2 + k] = 0;
        }
    }
    
    // Debug Sharp
    auto DebugSharp = [&]()
    {
        for (int i = 0; i < sharp_edges.size(); ++i) {
            if (sharp_edges[i] == 0)
                continue;
            int e = face_edgeIds[i/3][i%3];
            if (edge_diff[e][0] * edge_diff[e][1] != 0) {
                printf("Sharp condition violated!\n");
                exit(0);
            }
        }
        printf("Sharp condition pass!\n");
    };
    
    for (int i = 0; i < edge_diff.size(); ++i) {
        for (int j = 0; j < 2; ++j) {
            if (abs(edge_diff[i][j]) > 1) {
                edge_diff[i][j] /= abs(edge_diff[i][j]);
            }
        }
    }
    
#ifdef LOG_OUTPUT
    printf("Build Integer Constraints...\n");
#endif
    BuildIntegerConstraints();

    printf("MAX FLOW...\n");
    ComputeMaxFlow();

    // potential bug
#ifdef LOG_OUTPUT
    printf("subdivide...\n");
#endif
    subdivide_edgeDiff(F, V, N, Q, O, V2E, hierarchy.mE2E, boundary, nonManifold, edge_diff,
                       edge_values, face_edgeOrients, face_edgeIds, sharp_edges, singularities, 1);

    allow_changes.clear();
    allow_changes.resize(edge_diff.size() * 2, 1);
    for (int i = 0; i < sharp_edges.size(); ++i) {
        if (sharp_edges[i] == 0)
            continue;
        int e = face_edgeIds[i/3][i%3];
        for (int k = 0; k < 2; ++k) {
            if (edge_diff[e][k] == 0)
                allow_changes[e * 2 + k] = 0;
        }
    }

#ifdef LOG_OUTPUT
    printf("Fix flip advance...\n");
#endif

    DebugSharp();

    int t1 = GetCurrentTime64();
    FixFlipHierarchy();
    subdivide_edgeDiff(F, V, N, Q, O, V2E, hierarchy.mE2E, boundary, nonManifold, edge_diff,
                       edge_values, face_edgeOrients, face_edgeIds, sharp_edges, singularities, 1);
//    FixFlipSat();

    int t2 = GetCurrentTime64();
    printf("Flip use %lf\n", (t2 - t1) * 1e-3);

#ifdef LOG_OUTPUT
    printf("Post Linear Solver...\n");
#endif
    Optimizer::optimize_positions_fixed(hierarchy, edge_values, edge_diff, with_scale);
    AdvancedExtractQuad();
    FixValence();
    std::map<std::pair<int, int>, int> o2e;
    for (int i = 0; i < F_compact.size(); ++i) {
        for (int j = 0; j < 4; ++j) {
            int v1 = F_compact[i][j];
            int v2 = F_compact[i][(j + 1) % 4];
            o2e[std::make_pair(v1, v2)] = i * 4 + j;
        }
    }
    std::vector<std::vector<int> > v2o(V.cols());
    for (int i = 0; i < Vset.size(); ++i) {
        for (auto v : Vset[i]) {
            v2o[v].push_back(i);
        }
    }
    std::vector<Vector3d> diffs(F_compact.size() * 4, Vector3d(0, 0, 0));
    std::vector<int> diff_count(F_compact.size() * 4, 0);
    for (int i = 0; i < F.cols(); ++i) {
        for (int j = 0; j < 3; ++j) {
            int v1 = F(j, i);
            int v2 = F((j + 1) % 3, i);
            if (v1 != edge_values[face_edgeIds[i][j]].x)
                continue;
            if (edge_diff[face_edgeIds[i][j]].array().abs().sum() != 1)
                continue;
            if (v2o[v1].size() > 1 || v2o[v2].size() > 1)
                continue;
            for (auto o1 : v2o[v1]) {
                for (auto o2 : v2o[v2]) {
                    auto key = std::make_pair(o1, o2);
                    if (o2e.count(key)) {
                        int dedge = o2e[key];
                        Vector3d q_1 = Q.col(v1);
                        Vector3d q_2 = Q.col(v2);
                        Vector3d n_1 = N.col(v1);
                        Vector3d n_2 = N.col(v2);
                        Vector3d q_1_y = n_1.cross(q_1);
                        Vector3d q_2_y = n_2.cross(q_2);
                        auto index = compat_orientation_extrinsic_index_4(q_1, n_1, q_2, n_2);
                        //        double s_x1 = S(0, v1), s_y1 = S(1, v1);
                        //        double s_x2 = S(0, v2), s_y2 = S(1, v2);
                        int rank_diff = (index.second + 4 - index.first) % 4;
                        Vector3d qd_x = 0.5 * (rotate90_by(q_2, n_2, rank_diff) + q_1);
                        Vector3d qd_y = 0.5 * (rotate90_by(q_2_y, n_2, rank_diff) + q_1_y);
                        double scale_x = /*(with_scale ? 0.5 * (s_x1 + s_x2) : 1) */ hierarchy.mScale;
                        double scale_y = /*(with_scale ? 0.5 * (s_y1 + s_y2) : 1) */ hierarchy.mScale;
                        Vector2i diff = edge_diff[face_edgeIds[i][j]];
                        Vector3d C = diff[0] * scale_x * qd_x + diff[1] * scale_y * qd_y;
                        
                        diff_count[dedge] += 1;
                        diffs[dedge] += C;
                        auto key = std::make_pair(o2, o1);
                        if (o2e.count(key)) {
                            int dedge = o2e[key];
                            diff_count[dedge] += 1;
                            diffs[dedge] -= C;
                        }
                    }
                }
            }
        }
    }

    Optimizer::optimize_positions_dynamic(F, V, N, Q, Vset, O_compact, F_compact, V2E_compact,
                                          E2E_compact, sqrt(surface_area / F_compact.size()),
                                          diffs, diff_count, o2e);
    //    optimize_quad_positions(O_compact, N_compact, Q_compact, F_compact, V2E_compact,
    //    E2E_compact,
    //                            V, N, Q, O, F, V2E, hierarchy.mE2E, disajoint_tree,
    //                            hierarchy.mScale, false);
}
