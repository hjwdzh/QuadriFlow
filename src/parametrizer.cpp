#include "parametrizer.hpp"
#include "config.hpp"
#include "dedge.hpp"
#include "field-math.hpp"
#include "flow.hpp"
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
    ComputeOrientationSingularities();

    BuildEdgeInfo();

    for (int i = 0; i < face_edgeIds.size(); ++i) {
        for (int j = 0; j < 3; ++j) {
            if (face_edgeIds[i][j] == -1) {
                printf("OK, edge info is wrong!\n");
            }
        }
    }
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

    ComputeMaxFlow();

    printf("Analyze Valence...\n");
    // potential bug
#ifdef LOG_OUTPUT
    printf("subdivide...\n");
#endif
    subdivide_diff(F, V, N, Q, O, V2E, hierarchy.mE2E, boundary, nonManifold, edge_diff,
                   edge_values, face_edgeOrients, face_edgeIds, singularities, 1);

#ifdef LOG_OUTPUT
    printf("Fix flip advance...\n");
#endif
    int t1 = GetCurrentTime64();
    FixFlipHierarchy();

    int t2 = GetCurrentTime64();
    printf("Flip use %lf\n", (t2 - t1) * 1e-3);

//    subdivide_diff(F, V, N, Q, O, V2E, hierarchy.mE2E, boundary, nonManifold, edge_diff,
//                   edge_values, face_edgeOrients, face_edgeIds, singularities, 1);

    auto& E2E = hierarchy.mE2E;
    for (int i = 0; i < face_edgeIds.size(); ++i) {
        for (int j = 0; j < 3; ++j) {
            int e = i * 3 + j;
            int sum = 0;
            do {
                int e1 = e / 3 * 3 + (e + 2) % 3;
                e = E2E[e1];
                sum += abs(face_edgeOrients[e/3][e%3] + 6 - face_edgeOrients[e1/3][e1%3]);
            } while (e != i * 3 + j);
            if (sum % 4 != 0) {
                sing.insert(F(j, i));
            }
        }
    }
    std::queue<int> p;
    for (auto& s : sing)
        p.push(s);
    while (!p.empty()) {
        int e = V2E[p.front()];
        int e1 = e;
        do {
            Vector2i diff = edge_diff[face_edgeIds[e1/3][e1%3]];
            if (diff == Vector2i::Zero()) {
                int v = F((e1+1)%3, e1/3);
                if (sing.count(v) == 0) {
                    sing.insert(v);
                    p.push(v);
                }
            }
            e1 = e1 / 3 * 3 + (e1 + 2) % 3;
            e1 = E2E[e1];
        } while (e1 != e);
        p.pop();
    }
    int count = 0;
    for (int i = 0; i < face_edgeIds.size(); ++i) {
        Vector2i diff[3];
        for (int j = 0; j < 3; ++j) {
            diff[j] = rshift90(edge_diff[face_edgeIds[i][j]], face_edgeOrients[i][j]);
            if (abs(diff[j][0]) > 1 || abs(diff[j][1]) > 1) {
//                printf("Long edge... %d %d\n", diff[j][0], diff[j][1]);
            }
        }
        if (diff[0] + diff[1] + diff[2] != Vector2i::Zero()) {
            printf("Non zero!\n");
        }
        if (diff[0][0] * diff[1][1] - diff[0][1] * diff[1][0] < 0) {
            count += 1;
            flipped.insert(F(0, i));
            flipped.insert(F(1, i));
            flipped.insert(F(2, i));
        }
    }
    printf("Flipped %d\n", count);
    for (auto& s : flipped)
        p.push(s);
    while (!p.empty()) {
        int e = V2E[p.front()];
        int e1 = e;
        do {
            Vector2i diff = edge_diff[face_edgeIds[e1/3][e1%3]];
            if (diff == Vector2i::Zero()) {
                int v = F((e1+1)%3, e1/3);
                if (flipped.count(v) == 0) {
                    flipped.insert(v);
                    p.push(v);
                }
            }
            e1 = e1 / 3 * 3 + (e1 + 2) % 3;
            e1 = E2E[e1];
        } while (e1 != e);
        p.pop();
    }
    Optimizer::optimize_positions_fixed(hierarchy, edge_values, edge_diff, with_scale);

    ExtractQuadMesh();

#ifdef LOG_OUTPUT
    printf("Fix holes...\n");
#endif
    compute_direct_graph_quad(O_compact, F_compact, V2E_compact, E2E_compact, boundary_compact,
                              nonManifold_compact);
//    FixHoles();
    compute_direct_graph_quad(O_compact, F_compact, V2E_compact, E2E_compact, boundary_compact,
                              nonManifold_compact);
    // potential bug, not guarantee to have quads at holes!
#ifdef LOG_OUTPUT
    printf("Direct Quad Graph...\n");
#endif
#ifdef LOG_OUTPUT
    printf("Optimize quad positions...\n");
#endif

//    optimize_quad_positions(O_compact, N_compact, Q_compact, F_compact, V2E_compact, E2E_compact,
//                            V, N, Q, O, F, V2E, hierarchy.mE2E, disajoint_tree, hierarchy.mScale, false);
}
