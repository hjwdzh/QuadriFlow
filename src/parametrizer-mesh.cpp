#include "config.hpp"
#include "dedge.hpp"
#include "field-math.hpp"
#include "loader.hpp"
#include "merge-vertex.hpp"
#include "parametrizer.hpp"
#include "subdivide.hpp"

#include <queue>

void Parametrizer::Load(const char* filename) {
    load(filename, V, F);
    double maxV[3] = {-1e30, -1e30, -1e30};
    double minV[3] = {1e30, 1e30, 1e30};
#ifdef WITH_OMP
#pragma omp parallel for
#endif
    for (int i = 0; i < V.cols(); ++i) {
        for (int j = 0; j < 3; ++j) {
            maxV[j] = std::max(maxV[j], V(j, i));
            minV[j] = std::min(minV[j], V(j, i));
        }
    }
    double scale =
        std::max(std::max(maxV[0] - minV[0], maxV[1] - minV[1]), maxV[2] - minV[2]) * 0.5;
#ifdef WITH_OMP
#pragma omp parallel for
#endif
    for (int i = 0; i < V.cols(); ++i) {
        for (int j = 0; j < 3; ++j) {
            V(j, i) = (V(j, i) - (maxV[j] + minV[j]) * 0.5) / scale;
        }
    }
#ifdef LOG_OUTPUT
    printf("vertices size: %d\n", (int)V.cols());
    printf("faces size: %d\n", (int)F.cols());
#endif

    //    merge_close(V, F, 1e-6);
}

void Parametrizer::Initialize(int faces, int with_scale) {
    ComputeMeshStatus();
#ifdef PERFORMANCE_TEST
    num_vertices = V.cols() * 10;
    num_faces = num_vertices;
    scale = sqrt(surface_area / num_faces);
#else
    if (faces == -1) {
        num_vertices = V.cols();
        num_faces = num_vertices;
        scale = sqrt(surface_area / num_faces);
    } else {
        double face_area = surface_area / faces;
        num_vertices = faces;
        scale = std::sqrt(face_area) / 2;
    }
#endif
    double target_len = std::min(scale / 2, average_edge_length * 2);
#ifdef PERFORMANCE_TEST
    scale = sqrt(surface_area / V.cols());
#endif
    if (target_len < max_edge_length) {
        compute_direct_graph(V, F, V2E, E2E, boundary, nonManifold);
        subdivide(F, V, V2E, E2E, boundary, nonManifold, target_len);
    }
    int t1 = GetCurrentTime64();
    compute_direct_graph(V, F, V2E, E2E, boundary, nonManifold);
    generate_adjacency_matrix_uniform(F, V2E, E2E, nonManifold, adj);

    ComputeSmoothNormal();
    ComputeVertexArea();

    if (with_scale) {
        triangle_space.resize(F.cols());
#ifdef WITH_OMP
#pragma omp parallel for
#endif
        for (int i = 0; i < F.cols(); ++i) {
            Matrix3d p, q;
            p.col(0) = V.col(F(1, i)) - V.col(F(0, i));
            p.col(1) = V.col(F(2, i)) - V.col(F(0, i));
            p.col(2) = Nf.col(i);
            q = p.inverse();
            triangle_space[i].resize(2, 3);
            for (int j = 0; j < 2; ++j) {
                for (int k = 0; k < 3; ++k) {
                    triangle_space[i](j, k) = q(j, k);
                }
            }
        }
    }
#ifdef LOG_OUTPUT
    printf("V: %d F: %d\n", (int)V.cols(), (int)F.cols());
#endif
    hierarchy.mA[0] = std::move(A);
    hierarchy.mAdj[0] = std::move(adj);
    hierarchy.mN[0] = std::move(N);
    hierarchy.mV[0] = std::move(V);
    hierarchy.mE2E = std::move(E2E);
    hierarchy.mF = std::move(F);
    hierarchy.Initialize(scale, with_scale);
    int t2 = GetCurrentTime64();
    printf("Initialize use time: %lf\n", (t2 - t1) * 1e-3);
}

void Parametrizer::ComputeMeshStatus() {
    surface_area = 0;
    average_edge_length = 0;
    max_edge_length = 0;
    for (int f = 0; f < F.cols(); ++f) {
        Vector3d v[3] = {V.col(F(0, f)), V.col(F(1, f)), V.col(F(2, f))};
        double area = 0.5f * (v[1] - v[0]).cross(v[2] - v[0]).norm();
        surface_area += area;
        for (int i = 0; i < 3; ++i) {
            double len = (v[(i + 1) % 3] - v[i]).norm();
            average_edge_length += len;
            if (len > max_edge_length) max_edge_length = len;
        }
    }
    average_edge_length /= (F.cols() * 3);
}

void Parametrizer::ComputeSmoothNormal() {
    /* Compute face normals */
    Nf.resize(3, F.cols());
#ifdef WITH_OMP
#pragma omp parallel for
#endif
    for (int f = 0; f < F.cols(); ++f) {
        Vector3d v0 = V.col(F(0, f)), v1 = V.col(F(1, f)), v2 = V.col(F(2, f)),
                 n = (v1 - v0).cross(v2 - v0);
        double norm = n.norm();
        if (norm < RCPOVERFLOW) {
            n = Vector3d::UnitX();
        } else {
            n /= norm;
        }
        Nf.col(f) = n;
    }

    N.resize(3, V.cols());
#ifdef WITH_OMP
#pragma omp parallel for
#endif
    for (int i = 0; i < V2E.rows(); ++i) {
        int edge = V2E[i];
        if (nonManifold[i] || edge == -1) {
            N.col(i) = Vector3d::UnitX();
            continue;
        }

        int stop = edge;
        Vector3d normal = Vector3d::Zero();
        do {
            int idx = edge % 3;

            Vector3d d0 = V.col(F((idx + 1) % 3, edge / 3)) - V.col(i);
            Vector3d d1 = V.col(F((idx + 2) % 3, edge / 3)) - V.col(i);
            double angle = fast_acos(d0.dot(d1) / std::sqrt(d0.squaredNorm() * d1.squaredNorm()));

            /* "Computing Vertex Normals from Polygonal Facets"
             by Grit Thuermer and Charles A. Wuethrich, JGT 1998, Vol 3 */
            if (std::isfinite(angle)) normal += Nf.col(edge / 3) * angle;

            int opp = E2E[edge];
            if (opp == -1) break;

            edge = dedge_next_3(opp);
        } while (edge != stop);
        double norm = normal.norm();
        N.col(i) = norm > RCPOVERFLOW ? Vector3d(normal / norm) : Vector3d::UnitX();
    }
}

void Parametrizer::ComputeVertexArea() {
    A.resize(V.cols());
    A.setZero();

#ifdef WITH_OMP
#pragma omp parallel for
#endif
    for (int i = 0; i < V2E.size(); ++i) {
        int edge = V2E[i], stop = edge;
        if (nonManifold[i] || edge == -1) continue;
        double vertex_area = 0;
        do {
            int ep = dedge_prev_3(edge), en = dedge_next_3(edge);

            Vector3d v = V.col(F(edge % 3, edge / 3));
            Vector3d vn = V.col(F(en % 3, en / 3));
            Vector3d vp = V.col(F(ep % 3, ep / 3));

            Vector3d face_center = (v + vp + vn) * (1.0f / 3.0f);
            Vector3d prev = (v + vp) * 0.5f;
            Vector3d next = (v + vn) * 0.5f;

            vertex_area += 0.5f * ((v - prev).cross(v - face_center).norm() +
                                   (v - next).cross(v - face_center).norm());

            int opp = E2E[edge];
            if (opp == -1) break;
            edge = dedge_next_3(opp);
        } while (edge != stop);

        A[i] = vertex_area;
    }
}

void Parametrizer::ExtractQuadMesh() {
    auto& V = hierarchy.mV[0];
    auto& F = hierarchy.mF;
    auto& Q = hierarchy.mQ[0];
    auto& N = hierarchy.mN[0];
    auto& O = hierarchy.mO[0];

    disajoint_tree = DisajointTree(V.cols());
    for (int i = 0; i < edge_diff.size(); ++i) {
        if (edge_diff[i] == Vector2i::Zero()) {
            int vv0 = edge_values[i].x;
            int vv1 = edge_values[i].y;
            disajoint_tree.Merge(vv0, vv1);
        }
    }
    disajoint_tree.BuildCompactParent();

    int num_v = disajoint_tree.CompactNum();
    O_compact.resize(num_v, Vector3d::Zero());
    Q_compact.resize(num_v, Vector3d::Zero());
    N_compact.resize(num_v, Vector3d::Zero());
    counter.resize(num_v, 0);
    for (int i = 0; i < O.cols(); ++i) {
        int compact_v = disajoint_tree.Index(i);
        O_compact[compact_v] += O.col(i);
        N_compact[compact_v] = N_compact[compact_v] * counter[compact_v] + N.col(i);
        N_compact[compact_v].normalize();
        if (counter[compact_v] == 0)
            Q_compact[compact_v] = Q.col(i);
        else {
            auto pairs = compat_orientation_extrinsic_4(Q_compact[compact_v], N_compact[compact_v],
                                                        Q.col(i), N.col(i));
            Q_compact[compact_v] = (pairs.first * counter[compact_v] + pairs.second).normalized();
        }
        counter[compact_v] += 1;
    }
    for (int i = 0; i < O_compact.size(); ++i) {
        O_compact[i] /= counter[i];
    }

#ifdef LOG_OUTPUT
    printf("extract graph...\n");
#endif
    std::vector<std::set<int>> vertices(num_v), complete_set(num_v);
    for (int i = 0; i < edge_diff.size(); ++i) {
        int p1 = disajoint_tree.Index(edge_values[i].x);
        int p2 = disajoint_tree.Index(edge_values[i].y);
        if (p1 == p2) continue;
        complete_set[p1].insert(p2);
        complete_set[p2].insert(p1);
        if (abs(edge_diff[i][0]) + abs(edge_diff[i][1]) == 1) {
            vertices[p1].insert(p2);
            vertices[p2].insert(p1);
        }
    }

#ifdef LOG_OUTPUT
    printf("extract bad vertices...\n");
#endif
    bad_vertices.resize(num_v, 0);
    std::set<DEdge> bad_edges;
    
    std::queue<int> badq;
    for (int i = 0; i < num_v; ++i) {
        if (vertices[i].size() < 3) {
            badq.push(i);
            bad_vertices[i] = 1;
        }
    }
    while (!badq.empty()) {
        int v = badq.front();
        badq.pop();
        for (auto& v1 : vertices[v]) {
            vertices[v1].erase(v);
            if (vertices[v1].size() < 3 && bad_vertices[v1] == 0) {
                bad_vertices[v1] = 1;
                badq.push(v1);
            }
        }
    }
    /*
    for (int i = 0; i < F.cols(); ++i) {
        int v0 = F(0, i), p0 = disajoint_tree.Index(v0);
        int v1 = F(1, i), p1 = disajoint_tree.Index(v1);
        int v2 = F(2, i), p2 = disajoint_tree.Index(v2);
        if (p0 == p1 || p1 == p2 || p2 == p0) continue;
        Vector2i diff[3];
        for (int j = 0; j < 3; ++j) {
            int eid = face_edgeIds[i][j];
            diff[j] = rshift90(edge_diff[eid], face_edgeOrients[i][j]);
        }
        int a = -diff[0][0] * diff[2][1] + diff[0][1] * diff[2][0];
        if (a < 0) {
            for (int j = 0; j < 3; ++j) {
                int t1 = disajoint_tree.Index(F(j, i));
                int t2 = disajoint_tree.Index(F((j + 1) % 3, i));
                if (t1 != t2) bad_edges.insert(DEdge(t1, t2));
            }
        }
    }
    */
#ifdef LOG_OUTPUT
    printf("extract quad cells...\n");
#endif
    std::map<DEdge, std::pair<Vector3i, Vector3i>> quad_cells;
    for (int i = 0; i < F.cols(); ++i) {
        int v0 = F(0, i), p0 = disajoint_tree.Index(v0);
        int v1 = F(1, i), p1 = disajoint_tree.Index(v1);
        int v2 = F(2, i), p2 = disajoint_tree.Index(v2);
        if (p0 != p1 && p1 != p2 && p2 != p0 && bad_vertices[p0] == 0 && bad_vertices[p1] == 0 &&
            bad_vertices[p2] == 0 && bad_edges.count(DEdge(p0, p1)) == 0 &&
            bad_edges.count(DEdge(p1, p2)) == 0 && bad_edges.count(DEdge(p2, p0)) == 0) {
            auto diff1 = edge_diff[face_edgeIds[i][0]];
            auto diff2 = edge_diff[face_edgeIds[i][1]];
            auto diff3 = edge_diff[face_edgeIds[i][2]];
            int orient1 = face_edgeOrients[i][0];
            int orient2 = face_edgeOrients[i][2];
            auto d1 = rshift90(diff1, orient1);
            auto d2 = rshift90(-diff3, orient2);
//            if (d1[0] * d2[1] - d1[1] * d2[0] < 0) continue;
            DEdge eid;
            if (abs(diff1[0]) == 1 && abs(diff1[1]) == 1) {
                eid = DEdge(p0, p1);
            } else if (abs(diff2[0]) == 1 && abs(diff2[1]) == 1) {
                int t = p0;
                p0 = p1;
                p1 = p2;
                p2 = t;
                eid = DEdge(p0, p1);
            } else if (abs(diff3[0]) == 1 && abs(diff3[1]) == 1) {
                int t = p1;
                p1 = p0;
                p0 = p2;
                p2 = t;
                eid = DEdge(p0, p1);
            } else {
                continue;
            }
            if (quad_cells.count(eid) == 0)
                quad_cells[eid] = std::make_pair(Vector3i(p0, p1, p2), Vector3i(-100, -100, -100));
            else
                quad_cells[eid].second = Vector3i(p0, p1, p2);
        }
    }
#ifdef LOG_OUTPUT
    printf("extract quads...\n");
#endif
    for (auto& c : quad_cells) {
        if (c.second.second != Vector3i(-100, -100, -100)) {
            F_compact.push_back(Vector4i(c.second.first[0], c.second.second[2], c.second.first[1],
                                         c.second.first[2]));
        }
    }
}

void Parametrizer::OutputMesh(const char* obj_name) {
    std::vector<int> compact_answer(bad_vertices.size());
    compact_answer[0] = 1 - bad_vertices[0];
    for (int i = 1; i < bad_vertices.size(); ++i) {
        compact_answer[i] = compact_answer[i - 1] + (1 - bad_vertices[i]);
    }
    std::ofstream os(obj_name);
    for (int i = 0; i < bad_vertices.size(); ++i) {
        if (bad_vertices[i]) continue;
        os << "v " << O_compact[i][0] << " " << O_compact[i][1] << " " << O_compact[i][2] << "\n";
    }
    for (int i = 0; i < F_compact.size(); ++i) {
        os << "f " << compact_answer[F_compact[i][0]] << " " << compact_answer[F_compact[i][1]]
           << " " << compact_answer[F_compact[i][2]] << " " << compact_answer[F_compact[i][3]]
           << "\n";
    }
    os.close();
}
