#include "parametrizer.hpp"

#include "optimizer.hpp"
#include <vector>
#include <unordered_map>
#include <queue>

void Parametrizer::BuildEdgeInfo() {
    auto& F = hierarchy.mF;
    auto& E2E = hierarchy.mE2E;
    
    edge_diff.clear();
    edge_values.clear();
    face_edgeIds.resize(F.cols(), Vector3i(-1, -1, -1));
    for (int i = 0; i < F.cols(); ++i) {
        for (int j = 0; j < 3; ++j) {
            int k1 = j, k2 = (j + 1) % 3;
            int v1 = F(k1, i);
            int v2 = F(k2, i);
            DEdge e2(v1, v2);
            Vector2i diff2;
            int rank2;
            if (v1 > v2) {
                rank2 = pos_rank(k2, i);
                diff2 = rshift90(Vector2i(-pos_index(k1 * 2, i), -pos_index(k1 * 2 + 1, i)), rank2);
            } else {
                rank2 = pos_rank(k1, i);
                diff2 = rshift90(Vector2i(pos_index(k1 * 2, i), pos_index(k1 * 2 + 1, i)), rank2);
            }
            int current_eid = i * 3 + k1;
            int eid = E2E[current_eid];
            int eID1 = face_edgeIds[current_eid / 3][current_eid % 3];
            int eID2 = face_edgeIds[eid / 3][eid % 3];
            if (eID1 == -1) {
                eID2 = edge_values.size();
                edge_values.push_back(e2);
                edge_diff.push_back(diff2);
                face_edgeIds[i][k1] = eID2;
                face_edgeIds[eid / 3][eid % 3] = eID2;
            } else if (!singularities.count(i)) {
                edge_diff[eID2] = diff2;
            }
        }
    }
}

void Parametrizer::BuildIntegerConstraints() {
    auto& F = hierarchy.mF;
    auto& Q = hierarchy.mQ[0];
    auto& N = hierarchy.mN[0];
    face_edgeOrients.resize(F.cols());
    
    // undirected edge to direct edge
    std::vector<std::pair<int, int> > E2D(edge_diff.size(), std::make_pair(-1, -1));
    for (int i = 0; i < F.cols(); ++i) {
        int v0 = F(0, i);
        int v1 = F(1, i);
        int v2 = F(2, i);
        DEdge e0(v0, v1), e1(v1, v2), e2(v2, v0);
        const Vector3i& eid = face_edgeIds[i];
        Vector2i variable_id[3];
        for (int i = 0; i < 3; ++i) {
            variable_id[i] = Vector2i(eid[i] * 2 + 1, eid[i] * 2 + 2);
        }
        auto index1 = compat_orientation_extrinsic_index_4(Q.col(v0), N.col(v0), Q.col(v1), N.col(v1));
        auto index2 = compat_orientation_extrinsic_index_4(Q.col(v0), N.col(v0), Q.col(v2), N.col(v2));

        int rank1 = (index1.first - index1.second + 4) % 4; // v1 -> v0
        int rank2 = (index2.first - index2.second + 4) % 4; // v2 -> v0
        int orients[3] = {0}; // == {0, 0, 0}
        if (v1 < v0) {
            variable_id[0] = -rshift90(variable_id[0], rank1);
            orients[0] = (rank1 + 2) % 4;
        } else {
            orients[0] = 0;
        }
        if (v2 < v1) {
            variable_id[1] = -rshift90(variable_id[1], rank2);
            orients[1] = (rank2 + 2) % 4;
        } else {
            variable_id[1] = rshift90(variable_id[1], rank1);
            orients[1] = rank1;
        }
        if (v2 < v0) {
            variable_id[2] = rshift90(variable_id[2], rank2);
            orients[2] = rank2;
        } else {
            variable_id[2] = -variable_id[2];
            orients[2] = 2;
        }
        face_edgeOrients[i] = Vector3i(orients[0], orients[1], orients[2]);
        for (int j = 0; j < 3; ++j) {
            int eid = face_edgeIds[i][j];
            if (E2D[eid].first == -1)
                E2D[eid].first = i * 3 + j;
            else
                E2D[eid].second = i * 3 + j;
        }
    }

    // a face disajoint tree
    DisajointOrientTree disajoint_orient_tree = DisajointOrientTree(F.cols());
    // merge the whole face graph except for the singularity in which there exists a spanning tree which contains consistent orientation
    std::vector<int> sharpUE(E2D.size());
    for (int i = 0; i < sharp_edges.size(); ++i) {
        if (sharp_edges[i]) {
            sharpUE[face_edgeIds[i/3][i%3]] = 1;
        }
    }
    
    for (int i = 0; i < E2D.size(); ++i) {
        auto& edge_c = E2D[i];
        int f0 = edge_c.first / 3;
        int f1 = edge_c.second / 3;
        if (singularities.count(f0) || singularities.count(f1) || sharpUE[i]) continue;
        int orient1 = face_edgeOrients[f0][edge_c.first % 3];
        int orient0 = (face_edgeOrients[f1][edge_c.second % 3] + 2) % 4;
        disajoint_orient_tree.Merge(f0, f1, orient0, orient1);
    }
    
    for (int i = 0; i < sharpUE.size(); ++i) {
        if (sharpUE[i] == 0)
            continue;
        auto& edge_c = E2D[i];
        int f0 = edge_c.first / 3;
        int f1 = edge_c.second / 3;
        if (singularities.count(f0) || singularities.count(f1))
            continue;
        int orient1 = face_edgeOrients[f0][edge_c.first % 3];
        int orient0 = (face_edgeOrients[f1][edge_c.second % 3] + 2) % 4;
        disajoint_orient_tree.Merge(f0, f1, orient0, orient1);
    }
    
    // merge singularity later
    for (auto& f : singularities) {
        for (int i = 0; i < 3; ++i) {
            auto& edge_c = E2D[face_edgeIds[f.first][i]];
            int v0 = edge_c.first / 3;
            int v1 = edge_c.second / 3;
            int orient1 = face_edgeOrients[v0][edge_c.first % 3];
            int orient0 = (face_edgeOrients[v1][edge_c.second % 3] + 2) % 4;
            disajoint_orient_tree.Merge(v0, v1, orient0, orient1);
        }
    }
    
    // all the face has the same parent.  we rotate every face to the space of that parent.
    for (int i = 0; i < face_edgeOrients.size(); ++i) {
        for (int j = 0; j < 3; ++j) {
            face_edgeOrients[i][j] = (face_edgeOrients[i][j] + disajoint_orient_tree.Orient(i)) % 4;
        }
    }

    // each connected component are labeled with the same color
    std::vector<int> colors(face_edgeIds.size(), -1);
    int num_component = 0;
    // label the connected component connected by non-fixed edges
    // we need this because we need sink flow (demand) == source flow (supply) for each component rather than global
    for (int i = 0; i < colors.size(); ++i) {
        if (colors[i] != -1)
            continue;
        colors[i] = num_component;
        std::queue<int> q;
        q.push(i);
        int counter = 0;
        while (!q.empty()) {
            int v = q.front();
            q.pop();
            for (int i = 0; i < 3; ++i) {
                int e = face_edgeIds[v][i];
                if (abs(face_edgeOrients[E2D[e].first/3][E2D[e].first%3] -
                        face_edgeOrients[E2D[e].second/3][E2D[e].second%3] + 4) % 4 != 2) {
                    continue;
                }
                for (int k = 0; k < 2; ++k) {
                    int f = (k == 0) ? E2D[e].first/3 : E2D[e].second / 3;
                    if (colors[f] == -1) {
                        colors[f] = num_component;
                        q.push(f);
                    }
                }
            }
            counter += 1;
        }
        num_component += 1;
    }
    
    std::vector<int> total_flows(num_component);
    // check if each component is full-flow
    for (int i = 0; i < face_edgeIds.size(); ++i) {
        Vector2i diff(0, 0);
        for (int j = 0; j < 3; ++j) {
            int orient = face_edgeOrients[i][j];
            diff += rshift90(edge_diff[face_edgeIds[i][j]], orient);
        }
        total_flows[colors[i]] += diff[0] + diff[1];
    }

    // build "variable"
    variables.resize(edge_diff.size() * 2, std::make_pair(Vector2i(-1, -1), 0));
    for (int i = 0; i < face_edgeIds.size(); ++i) {
        for (int j = 0; j < 3; ++j) {
            Vector2i sign = rshift90(Vector2i(1, 1), face_edgeOrients[i][j]);
            int eid = face_edgeIds[i][j];
            Vector2i index = rshift90(Vector2i(eid * 2, eid * 2 + 1), face_edgeOrients[i][j]);
            for (int k = 0; k < 2; ++k) {
                auto& p = variables[abs(index[k])];
                if (p.first[0] == -1)
                    p.first[0] = i * 2 + k;
                else
                    p.first[1] = i * 2 + k;
                p.second += sign[k];
            }
        }
    }

    // fixed variable that might be manually modified.
    // modified_variables[component_od][].first = fixed_variable_id
    // modified_variables[component_od][].second = 1 if two positive signs -1 if two negative signs
    std::vector<std::vector<std::pair<int, int> > > modified_variables(total_flows.size());
    for (int i = 0; i < variables.size(); ++i) {
        if (variables[i].second != 0 && allow_changes[i] == 1) {
            int find = colors[variables[i].first[0]/2];
            if (total_flows[find] > 0) {
                if (variables[i].second > 0 && edge_diff[i / 2][i % 2] > -1) {
                    modified_variables[find].push_back(std::make_pair(i, -1));
                }
                if (variables[i].second < 0 && edge_diff[i / 2][i % 2] < 1) {
                    modified_variables[find].push_back(std::make_pair(i, 1));
                }
            } else if (total_flows[find] < 0) {
                if (variables[i].second < 0 && edge_diff[i / 2][i % 2] > -1) {
                    modified_variables[find].push_back(std::make_pair(i, -1));
                }
                if (variables[i].second > 0 && edge_diff[i / 2][i % 2] < 1) {
                    modified_variables[find].push_back(std::make_pair(i, 1));
                }
            }
        }
    }

    // uniformly random manually modify variables so that the network has full flow.
    for (auto& modified_var : modified_variables)
        std::random_shuffle(modified_var.begin(), modified_var.end());

    for (int j = 0; j < total_flows.size(); ++j) {
        if (total_flows[j] == 0)
            continue;
        for (int i = 0; i < abs(total_flows[j]) / 2; ++i) {
            auto& info = modified_variables[j][i];
            edge_diff[info.first / 2][info.first % 2] += info.second;
        }
    }
    
    total_flows.clear();
    total_flows.resize(num_component, 0);
    // check if each component is full-flow
    for (int i = 0; i < face_edgeIds.size(); ++i) {
        Vector2i diff(0, 0);
        for (int j = 0; j < 3; ++j) {
            int orient = face_edgeOrients[i][j];
            diff += rshift90(edge_diff[face_edgeIds[i][j]], orient);
        }
        total_flows[colors[i]] += diff[0] + diff[1];
    }
    
    printf("Finish...\n");
}

void Parametrizer::ComputeMaxFlow() {
    printf("Downsample!\n");
    hierarchy.DownsampleEdgeGraph(face_edgeOrients, face_edgeIds, edge_diff, sharp_edges, allow_changes, 6);
    for (int i = 0; i < hierarchy.mEdgeDiff.size(); ++i) {
        printf("sample %d %d\n", hierarchy.mEdgeDiff[i].size(), hierarchy.mAllowChanges[i].size());
    }
    Optimizer::optimize_integer_constraints(hierarchy, singularities);
    hierarchy.UpdateGraphValue(face_edgeOrients, face_edgeIds, edge_diff);
}
