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
    face_edgeOrients.resize(F.cols());
    
    std::vector<std::pair<int, int> > E2F(edge_diff.size(), std::make_pair(-1, -1));
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
        for (int j = 0; j < 3; ++j) {
            int eid = face_edgeIds[i][j];
            if (E2F[eid].first == -1)
                E2F[eid].first = i * 3 + j;
            else
                E2F[eid].second = i * 3 + j;
        }
    }
    
    DisajointOrientTree disajoint_orient_tree = DisajointOrientTree(F.cols());
    for (int i = 0; i < E2F.size(); ++i) {
        auto& edge_c = E2F[i];
        int v0 = edge_c.first / 3;
        int v1 = edge_c.second / 3;
        if (singularities.count(v0) || singularities.count(v1)) continue;
        int orient1 = face_edgeOrients[v0][edge_c.first % 3];
        int orient0 = (face_edgeOrients[v1][edge_c.second % 3] + 2) % 4;
        disajoint_orient_tree.Merge(v0, v1, orient0, orient1);
    }
    
    for (auto& f : singularities) {
        for (int i = 0; i < 3; ++i) {
            auto& edge_c = E2F[face_edgeIds[f.first][i]];
            int v0 = edge_c.first / 3;
            int v1 = edge_c.second / 3;
            int orient1 = face_edgeOrients[v0][edge_c.first % 3];
            int orient0 = (face_edgeOrients[v1][edge_c.second % 3] + 2) % 4;
            disajoint_orient_tree.Merge(v0, v1, orient0, orient1);
        }
    }
    
    for (int i = 0; i < face_edgeOrients.size(); ++i) {
        for (int j = 0; j < 3; ++j) {
            face_edgeOrients[i][j] = (face_edgeOrients[i][j] + disajoint_orient_tree.Orient(i)) % 4;
        }
    }
    
    for (auto& f : singularities) {
        Vector3i cuts[4];
        for (int j = 0; j < 4; ++j) {
            for (int i = 0; i < 3; ++i) {
                int e = face_edgeIds[f.first][i];
                int orient1 = face_edgeOrients[E2F[e].first/3][E2F[e].first%3];
                int orient2 = face_edgeOrients[E2F[e].second/3][E2F[e].second%3];
                if (E2F[e].first / 3 == f.first)
                    orient1 += j;
                else
                    orient2 += j;
                if (abs(orient1 - orient2 + 10) % 4 != 0)
                    cuts[j][i] = 1;
                else
                    cuts[j][i] = 0;
            }
        }
        int min_cuts = 3;
        int ind = 4;
        for (int j = 0; j < 4; ++j) {
            int cut = cuts[j][0] + cuts[j][1] + cuts[j][2];
            if (cut < min_cuts) {
                min_cuts = cut;
                ind = j;
            }
        }
        int res = 0;
        if (cuts[ind][0] != 0 || cuts[ind][1] != 0 || cuts[ind][2] != 0) {
            while (cuts[ind][res] == 0)
                res += 1;
        }
        for (int j = 0; j < 3; ++j) {
            if (j == res)
                face_edgeOrients[f.first][j] = (face_edgeOrients[f.first][j] + f.second + ind) % 4;
            else
                face_edgeOrients[f.first][j] = (face_edgeOrients[f.first][j] + ind) % 4;
        }
    }
    
    std::vector<int> colors(face_edgeIds.size(), -1);
    int num_v = 0;
    for (int i = 0; i < colors.size(); ++i) {
        if (colors[i] != -1)
            continue;
        colors[i] = num_v;
        std::queue<int> q;
        q.push(i);
        int counter = 0;
        while (!q.empty()) {
            int v = q.front();
            q.pop();
            for (int i = 0; i < 3; ++i) {
                int e = face_edgeIds[v][i];
                    if (abs(face_edgeOrients[E2F[e].first/3][E2F[e].first%3] -
                            face_edgeOrients[E2F[e].second/3][E2F[e].second%3] + 4) % 4 != 2) {
                        continue;
                    }
                for (int k = 0; k < 2; ++k) {
                    int f = (k == 0) ? E2F[e].first/3 : E2F[e].second / 3;
                    if (colors[f] == -1) {
                        colors[f] = num_v;
                        q.push(f);
                    }
                }
            }
            counter += 1;
        }
        num_v += 1;
    }
    
    DisajointTree segments(num_v);
    for (auto& f : singularities) {
        for (int i = 0; i < 3; ++i) {
            int eid = face_edgeIds[f.first][i];
            int f1 = E2F[eid].first/3;
            int f2 = E2F[eid].second/3;
            if (segments.Parent(colors[f1]) != segments.Parent(colors[f2])) {
                segments.Merge(colors[f1], colors[f2]);
                if (f1 == f.first) {
                    face_edgeOrients[f1][E2F[eid].first%3] = (face_edgeOrients[f2][E2F[eid].second%3] + 2) % 4;
                } else {
                    face_edgeOrients[f2][E2F[eid].second%3] = (face_edgeOrients[f1][E2F[eid].first%3] + 2) % 4;
                }
            }
        }
    }
    
    for (int eid = 0; eid < E2F.size(); ++eid) {
        int f1 = E2F[eid].first/3;
        int f2 = E2F[eid].second/3;
        if (segments.Parent(colors[f1]) != segments.Parent(colors[f2])) {
            segments.Merge(colors[f1], colors[f2]);
            face_edgeOrients[f1][E2F[eid].first%3] = (face_edgeOrients[f2][E2F[eid].second%3] + 2) % 4;
        }
    }
    
    segments.BuildCompactParent();
    
    std::vector<int> total_flows(segments.CompactNum());
    for (int i = 0; i < face_edgeIds.size(); ++i) {
        Vector2i diff(0, 0);
        for (int j = 0; j < 3; ++j) {
            int orient = face_edgeOrients[i][j];
            diff += rshift90(edge_diff[face_edgeIds[i][j]], orient);
        }
        total_flows[segments.Index(colors[i])] += diff[0] + diff[1];
    }

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
    cuts.clear();

    std::vector<std::vector<std::pair<int, int> > > modified_variables(total_flows.size());
    for (int i = 0; i < variables.size(); ++i) {
        if (variables[i].second != 0) {
            int find = segments.Index(colors[variables[i].first[0]/2]);
            cuts.insert(edge_values[i / 2]);
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

    for (auto& modified_var : modified_variables)
        std::random_shuffle(modified_var.begin(), modified_var.end());
    for (int j = 0; j < total_flows.size(); ++j) {
        for (int i = 0; i < abs(total_flows[j]) / 2; ++i) {
            auto& info = modified_variables[j][i];
            edge_diff[info.first / 2][info.first % 2] += info.second;
        }
    }
    
    for (int i = 0; i < segments.CompactNum(); ++i) {
        int sum = 0;
        for (int j = 0; j < face_edgeIds.size(); ++j) {
            if (segments.Index(colors[j]) != i)
                continue;
            Vector2i diff(0, 0);
            for (int k = 0; k < 3; ++k) {
                diff += rshift90(edge_diff[face_edgeIds[j][k]], face_edgeOrients[j][k]);
            }
            sum += diff[0] + diff[1];
        }
    }
}

void Parametrizer::ComputeMaxFlow() {
    hierarchy.DownsampleEdgeGraph(face_edgeOrients, face_edgeIds, edge_diff, 6);
    Optimizer::optimize_integer_constraints(hierarchy, singularities);
    hierarchy.UpdateGraphValue(face_edgeOrients, face_edgeIds, edge_diff);
}

