#include "parametrizer.hpp"

#include "optimizer.hpp"
#include <vector>
#include <unordered_map>

void Parametrizer::BuildEdgeInfo() {
    auto& F = hierarchy.mF;
    auto& E2E = hierarchy.mE2E;
    
    edge_diff.clear();
    edge_values.clear();
    face_edgeIds.resize(F.cols(), Vector3i(-1, -1, -1));
    for (int i = 0; i < F.cols(); ++i) {
        for (int j = 0; j < 3; ++j) {
            //            if (face_edgeIds[i][j] != -1)
            //                continue;
            int k1 = j, k2 = (j + 1) % 3;
            int v1 = F(k1, i);
            int v2 = F(k2, i);
            DEdge e2(v1, v2);
            Vector2i diff2;
            int rank2;
            if (v1 > v2) {
                rank2 = pos_rank(k2, i);
                diff2 =
                rshift90(Vector2i(-pos_index(k1 * 2, i), -pos_index(k1 * 2 + 1, i)), rank2);
            } else {
                rank2 = pos_rank(k1, i);
                diff2 = rshift90(Vector2i(pos_index(k1 * 2, i), pos_index(k1 * 2 + 1, i)), rank2);
            }
            int current_eid = i * 3 + k1;
            int eid = E2E[current_eid];
            int eID2 = face_edgeIds[eid / 3][eid % 3];
            if (eID2 == -1) {
                eID2 = edge_values.size();
                edge_values.push_back(e2);
                edge_diff.push_back(diff2);
                face_edgeIds[i][k1] = eID2;
                if (eid != -1) face_edgeIds[eid / 3][eid % 3] = eID2;
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
    std::vector<Vector2i> sign_indices;
    face_edgeOrients.resize(F.cols());
    std::vector<Vector4i> edge_to_constraints;
    edge_to_constraints.resize(edge_values.size(), Vector4i(-1, -1, -1, -1));
    
    for (int i = 0; i < F.cols(); ++i) {
        int v0 = F(0, i);
        int v1 = F(1, i);
        int v2 = F(2, i);
        DEdge e0(v0, v1), e1(v1, v2), e2(v2, v0);
        const Vector3i& eid = face_edgeIds[i];
        Vector2i vid[3];
        for (int i = 0; i < 3; ++i) {
            vid[i] = Vector2i(eid[i] * 2 + 1, eid[i] * 2 + 2);
        }
        auto index1 =
        compat_orientation_extrinsic_index_4(Q.col(v0), N.col(v0), Q.col(v1), N.col(v1));
        auto index2 =
        compat_orientation_extrinsic_index_4(Q.col(v0), N.col(v0), Q.col(v2), N.col(v2));
        int rank1 = (index1.first - index1.second + 4) % 4;
        int rank2 = (index2.first - index2.second + 4) % 4;
        int orients[3] = {0};
        if (v1 < v0) {
            vid[0] = -rshift90(vid[0], rank1);
            orients[0] = (rank1 + 2) % 4;
        }
        if (v2 < v1) {
            vid[1] = -rshift90(vid[1], rank2);
            orients[1] = (rank2 + 2) % 4;
        } else {
            vid[1] = rshift90(vid[1], rank1);
            orients[1] = rank1;
        }
        if (v2 < v0) {
            vid[2] = rshift90(vid[2], rank2);
            orients[2] = rank2;
        } else {
            vid[2] = -vid[2];
            orients[2] = 2;
        }
        face_edgeOrients[i] = Vector3i(orients[0], orients[1], orients[2]);
        edge_to_constraints[eid[0]][(v0 > v1) * 2] = i;
        edge_to_constraints[eid[0]][(v0 > v1) * 2 + 1] = orients[0];
        edge_to_constraints[eid[1]][(v1 > v2) * 2] = i;
        edge_to_constraints[eid[1]][(v1 > v2) * 2 + 1] = orients[1];
        edge_to_constraints[eid[2]][(v2 > v0) * 2] = i;
        edge_to_constraints[eid[2]][(v2 > v0) * 2 + 1] = orients[2];
        
        for (int k = 0; k < 3; ++k) {
            sign_indices.push_back(vid[k]);
        }
    }
    
    DisajointOrientTree disajoint_orient_tree = DisajointOrientTree(F.cols());
    for (int i = 0; i < edge_to_constraints.size(); ++i) {
        auto& edge_c = edge_to_constraints[i];
        int v0 = edge_c[0];
        int v1 = edge_c[2];
        if (singularities.count(v0) || singularities.count(v1)) continue;
        int orient1 = edge_c[1];
        int orient0 = (edge_c[3] + 2) % 4;
        disajoint_orient_tree.Merge(v0, v1, orient0, orient1);
    }
    
    std::vector<Vector3i> sing_diff;
    std::vector<std::unordered_map<int, std::pair<int, int>>> sing_maps;
    std::vector<Vector3i> sing_orients;
    for (int i = 0; i < sign_indices.size(); i += 3) {
        int f = i / 3;
        int orient = disajoint_orient_tree.Orient(f);
        for (int j = 0; j < 3; ++j) {
            sign_indices[i + j] = rshift90(sign_indices[i + j], orient);
        }
        for (int j = 0; j < 2; ++j) {
            Vector3i sign, ind;
            for (int k = 0; k < 3; ++k) {
                ind[k] = abs(sign_indices[i + k][j]);
                if (ind[k] == 0) {
                    printf("OMG!\n");
                    exit(0);
                }
                sign[k] = sign_indices[i + k][j] / ind[k];
                ind[k] -= 1;
            }
            constraints_index.push_back(ind);
            constraints_sign.push_back(sign);
        }
        if (singularities.count(f)) {
            int orient_base = singularities[f];
            Vector3i diffs;
            Vector3i orient_diffs;
            for (int j = 0; j < 3; ++j) {
                int eid = face_edgeIds[f][(j + 1) % 3];
                int v0 = edge_to_constraints[eid][0];
                int v1 = edge_to_constraints[eid][2];
                int orientp0 = disajoint_orient_tree.Orient(v0) + edge_to_constraints[eid][1];
                int orientp1 = disajoint_orient_tree.Orient(v1) + edge_to_constraints[eid][3];
                int orient_diff = 0;
                if (v1 == f)
                    orient_diff = (orientp0 - orientp1 + 6) % 4;
                else
                    orient_diff = (orientp1 - orientp0 + 6) % 4;
                Vector2i sign_index[3];
                sign_index[0] = rshift90(sign_indices[i + j], (orient_base + orient_diff) % 4);
                sign_index[1] = rshift90(sign_indices[i + (j + 1) % 3], orient_diff);
                sign_index[2] = rshift90(sign_indices[i + (j + 2) % 3], orient_diff);
                int total_diff = 0;
                for (int k = 0; k < 2; ++k) {
                    auto ind = constraints_index[f * 2 + k];
                    auto sign = constraints_sign[f * 2 + k];
                    for (int l = 0; l < 3; ++l) {
                        ind[l] = abs(sign_index[l][k]);
                        sign[l] = sign_index[l][k] / ind[l];
                        ind[l] -= 1;
                    }
                    int diff1 = edge_diff[ind[0] / 2][ind[0] % 2];
                    int diff2 = edge_diff[ind[1] / 2][ind[1] % 2];
                    int diff3 = edge_diff[ind[2] / 2][ind[2] % 2];
                    int diff = sign[0] * diff1 + sign[1] * diff2 + sign[2] * diff3;
                    total_diff += diff;
                }
                orient_diffs[j] = orient_diff;
                diffs[j] = total_diff;
            }
            sing_diff.push_back(diffs);
            sing_orients.push_back(orient_diffs);
        }
    }
    
    int total_flow = 0;
    for (int i = 0; i < constraints_index.size(); ++i) {
        if (singularities.count(i / 2)) {
            continue;
        }
        auto index = constraints_index[i];
        auto sign = constraints_sign[i];
        int diff1 = edge_diff[index[0] / 2][index[0] % 2];
        int diff2 = edge_diff[index[1] / 2][index[1] % 2];
        int diff3 = edge_diff[index[2] / 2][index[2] % 2];
        int diff = sign[0] * diff1 + sign[1] * diff2 + sign[2] * diff3;
        total_flow += diff;
    }
    
    sing_maps.resize(sing_diff.size() + 1);
    sing_maps[0][total_flow] = std::make_pair(0, 0);
    for (int i = 0; i < sing_diff.size(); ++i) {
        auto& prev = sing_maps[i];
        auto& next = sing_maps[i + 1];
        for (auto& p : prev) {
            for (int j = 0; j < 3; ++j) {
                int v = p.first + sing_diff[i][j];
                int t = p.second.first + abs(sing_diff[i][j]);
                auto it = next.find(v);
                if (it == next.end())
                    next[v] = std::make_pair(t, j);
                else if (t < it->second.first)
                    it->second = std::make_pair(t, j);
            }
        }
    }
    
    std::vector<int> sing_selection;
    int target_flow = 0;
    while (sing_maps.back().count(target_flow) == 0 && sing_maps.back().count(-target_flow) == 0) {
        target_flow += 2;
    }
    if (sing_maps.back().count(target_flow) == 0) target_flow = -target_flow;
    int remain_flow = target_flow;
    for (int i = sing_diff.size(); i > 0; i--) {
        auto p = sing_maps[i][remain_flow];
        remain_flow -= sing_diff[i - 1][p.second];
        sing_selection.push_back(p.second);
    }
    std::reverse(sing_selection.begin(), sing_selection.end());
    int sing_count = 0;
    for (auto& f : singularities) {
        int select = sing_selection[sing_count];
        int orient_diff = sing_orients[sing_count++][select];
        auto& index1 = constraints_index[f.first * 2];
        auto& index2 = constraints_index[f.first * 2 + 1];
        auto& sign1 = constraints_sign[f.first * 2];
        auto& sign2 = constraints_sign[f.first * 2 + 1];
        
        int eid0;
        for (int i = 0; i < 3; ++i) {
            auto diff = Vector2i(sign1[i] * (index1[i] + 1), sign2[i] * (index2[i] + 1));
            int t = orient_diff;
            if (i == select) t = (t + f.second) % 4;
            int v0 = F(i, f.first);
            int v1 = F((i + 1) % 3, f.first);
            int eid = face_edgeIds[f.first][i];
            if ((select + 1) % 3 == i) eid0 = eid;
            edge_to_constraints[eid][(v0 > v1) * 2] = f.first;
            edge_to_constraints[eid][(v0 > v1) * 2 + 1] =
            (edge_to_constraints[eid][(v0 > v1) * 2 + 1] + t) % 4;
            face_edgeOrients[f.first][i] = (face_edgeOrients[f.first][i] + t) % 4;
            
            diff = rshift90(diff, t);
            index1[i] = abs(diff[0]);
            sign1[i] = diff[0] / abs(diff[0]);
            index1[i] -= 1;
            index2[i] = abs(diff[1]);
            sign2[i] = diff[1] / abs(diff[1]);
            index2[i] -= 1;
        }
        auto& edge_c = edge_to_constraints[eid0];
        int v0 = edge_c[0];
        int v1 = edge_c[2];
        
        int orient1 = edge_c[1];
        int orient0 = (edge_c[3] + 2) % 4;
        
        disajoint_orient_tree.Merge(v0, v1, orient0, orient1);
    }
    total_flow = 0;
    
    variables.resize(edge_diff.size() * 2, std::make_pair(Vector2i(-1, -1), 0));
    for (int i = 0; i < constraints_index.size(); ++i) {
        auto index = constraints_index[i];
        auto sign = constraints_sign[i];
        int diff1 = edge_diff[index[0] / 2][index[0] % 2];
        int diff2 = edge_diff[index[1] / 2][index[1] % 2];
        int diff3 = edge_diff[index[2] / 2][index[2] % 2];
        int diff = sign[0] * diff1 + sign[1] * diff2 + sign[2] * diff3;
        total_flow += diff;
        for (int j = 0; j < 3; ++j) {
            auto& p = variables[index[j]].first;
            if (sign[j] > 0)
                p[0] = i;
            else
                p[1] = i;
            variables[index[j]].second += sign[j];
        }
    }
    cuts.clear();
    std::vector<std::pair<int, int>> modified_variables;
    for (int i = 0; i < variables.size(); ++i) {
        if (variables[i].second != 0) {
            cuts.insert(edge_values[i / 2]);
            if (target_flow > 0) {
                if (variables[i].second > 0 && edge_diff[i / 2][i % 2] > -1) {
                    modified_variables.push_back(std::make_pair(i, -1));
                }
                if (variables[i].second < 0 && edge_diff[i / 2][i % 2] < 1) {
                    modified_variables.push_back(std::make_pair(i, 1));
                }
            } else if (target_flow < 0) {
                if (variables[i].second < 0 && edge_diff[i / 2][i % 2] > -1) {
                    modified_variables.push_back(std::make_pair(i, -1));
                }
                if (variables[i].second > 0 && edge_diff[i / 2][i % 2] < 1) {
                    modified_variables.push_back(std::make_pair(i, 1));
                }
            }
        }
    }
    
    std::random_shuffle(modified_variables.begin(), modified_variables.end());
    
    for (int i = 0; i < abs(target_flow) / 2; ++i) {
        auto& info = modified_variables[i];
        edge_diff[info.first / 2][info.first % 2] += info.second;
    }
    
    for (int i = 0; i < face_edgeOrients.size(); ++i) {
        for (int j = 0; j < 3; ++j) {
            face_edgeOrients[i][j] =
            (face_edgeOrients[i][j] + disajoint_orient_tree.Orient(i)) % 4;
        }
    }
}

void Parametrizer::ComputeMaxFlow() {
    hierarchy.DownsampleEdgeGraph(face_edgeOrients, face_edgeIds, edge_diff, 6);
    Optimizer::optimize_integer_constraints(hierarchy, singularities);
    hierarchy.UpdateGraphValue(face_edgeOrients, face_edgeIds, edge_diff);
}
