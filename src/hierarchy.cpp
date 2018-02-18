#include "hierarchy.hpp"
#include <fstream>
#include <unordered_map>
#include "config.hpp"
#include "field-math.hpp"
#ifdef WITH_TBB
#include "tbb_common.h"
#endif
#include <queue>
#include "pcg32/pcg32.h"
#include "pss/parallel_stable_sort.h"
Hierarchy::Hierarchy() {
    mAdj.resize(MAX_DEPTH + 1);
    mV.resize(MAX_DEPTH + 1);
    mN.resize(MAX_DEPTH + 1);
    mA.resize(MAX_DEPTH + 1);
    mPhases.resize(MAX_DEPTH + 1);
    mToLower.resize(MAX_DEPTH);
    mToUpper.resize(MAX_DEPTH);
}

#undef max

void Hierarchy::Initialize(double scale, int with_scale) {
    this->with_scale = with_scale;
    generate_graph_coloring_deterministic(mAdj[0], mV[0].cols(), mPhases[0]);
    for (int i = 0; i < MAX_DEPTH; ++i) {
        DownsampleGraph(mAdj[i], mV[i], mN[i], mA[i], mV[i + 1], mN[i + 1], mA[i + 1], mToUpper[i],
                        mToLower[i], mAdj[i + 1]);
        generate_graph_coloring_deterministic(mAdj[i + 1], mV[i + 1].cols(), mPhases[i + 1]);
        if (mV[i + 1].cols() == 1) {
            mAdj.resize(i + 2);
            mV.resize(i + 2);
            mN.resize(i + 2);
            mA.resize(i + 2);
            mToUpper.resize(i + 1);
            mToLower.resize(i + 1);
            break;
        }
    }
    mQ.resize(mV.size());
    mO.resize(mV.size());
    
    if (with_scale) {
        mS.resize(mV.size());
        mK.resize(mV.size());
    }
    
    mScale = scale;
#ifdef WITH_OMP
#pragma omp parallel for
#endif
    for (int i = 0; i < mV.size(); ++i) {
        mQ[i].resize(mN[i].rows(), mN[i].cols());
        mO[i].resize(mN[i].rows(), mN[i].cols());
        if (with_scale) {
            mS[i].resize(2, mN[i].cols());
            mK[i].resize(2, mN[i].cols());
        }
        for (int j = 0; j < mN[i].cols(); ++j) {
            Vector3d s, t;
            coordinate_system(mN[i].col(j), s, t);
            double angle = ((double)rand()) / RAND_MAX * 2 * M_PI;
            double x = ((double)rand()) / RAND_MAX * 2 - 1.f;
            double y = ((double)rand()) / RAND_MAX * 2 - 1.f;
            mQ[i].col(j) = s * std::cos(angle) + t * std::sin(angle);
            mO[i].col(j) = mV[i].col(j) + (s * x + t * y) * scale;
            if (with_scale) {
                mS[i].col(j) = Vector2d(1.0f, 1.0f);
                mK[i].col(j) = Vector2d(0.0, 0.0);
            }
        }
    }
#ifdef WITH_CUDA
    printf("copy to device...\n");
    CopyToDevice();
    printf("copy to device finish...\n");
#endif
}

#ifdef WITH_TBB
void Hierarchy::generate_graph_coloring_deterministic(const AdjacentMatrix& adj, int size,
                                                      std::vector<std::vector<int>>& phases) {
    struct ColorData {
        uint8_t nColors;
        uint32_t nNodes[256];
        ColorData() : nColors(0) {}
    };
    
    const uint8_t INVALID_COLOR = 0xFF;
    phases.clear();
    
    /* Generate a permutation */
    std::vector<uint32_t> perm(size);
    std::vector<tbb::spin_mutex> mutex(size);
    for (uint32_t i = 0; i < size; ++i) perm[i] = i;
    
    tbb::parallel_for(tbb::blocked_range<uint32_t>(0u, size, GRAIN_SIZE),
                      [&](const tbb::blocked_range<uint32_t>& range) {
                          pcg32 rng;
                          rng.advance(range.begin());
                          for (uint32_t i = range.begin(); i != range.end(); ++i) {
                              uint32_t j = i, k = rng.nextUInt(size - i) + i;
                              if (j == k) continue;
                              if (j > k) std::swap(j, k);
                              tbb::spin_mutex::scoped_lock l0(mutex[j]);
                              tbb::spin_mutex::scoped_lock l1(mutex[k]);
                              std::swap(perm[j], perm[k]);
                          }
                      });
    
    std::vector<uint8_t> color(size, INVALID_COLOR);
    ColorData colorData = tbb::parallel_reduce(
                                               tbb::blocked_range<uint32_t>(0u, size, GRAIN_SIZE), ColorData(),
                                               [&](const tbb::blocked_range<uint32_t>& range, ColorData colorData) -> ColorData {
                                                   std::vector<uint32_t> neighborhood;
                                                   bool possible_colors[256];
                                                   
                                                   for (uint32_t pidx = range.begin(); pidx != range.end(); ++pidx) {
                                                       uint32_t i = perm[pidx];
                                                       
                                                       neighborhood.clear();
                                                       neighborhood.push_back(i);
                                                       //            for (const Link *link = adj[i]; link != adj[i + 1]; ++link)
                                                       for (auto& link : adj[i]) neighborhood.push_back(link.id);
                                                       std::sort(neighborhood.begin(), neighborhood.end());
                                                       for (uint32_t j : neighborhood) mutex[j].lock();
                                                       
                                                       std::fill(possible_colors, possible_colors + colorData.nColors, true);
                                                       
                                                       //            for (const Link *link = adj[i]; link != adj[i + 1]; ++link) {
                                                       for (auto& link : adj[i]) {
                                                           uint8_t c = color[link.id];
                                                           if (c != INVALID_COLOR) {
                                                               while (c >= colorData.nColors) {
                                                                   possible_colors[colorData.nColors] = true;
                                                                   colorData.nNodes[colorData.nColors] = 0;
                                                                   colorData.nColors++;
                                                               }
                                                               possible_colors[c] = false;
                                                           }
                                                       }
                                                       
                                                       uint8_t chosen_color = INVALID_COLOR;
                                                       for (uint8_t j = 0; j < colorData.nColors; ++j) {
                                                           if (possible_colors[j]) {
                                                               chosen_color = j;
                                                               break;
                                                           }
                                                       }
                                                       if (chosen_color == INVALID_COLOR) {
                                                           if (colorData.nColors == INVALID_COLOR - 1)
                                                               throw std::runtime_error(
                                                                                        "Ran out of colors during graph coloring! "
                                                                                        "The input mesh is very likely corrupt.");
                                                           colorData.nNodes[colorData.nColors] = 1;
                                                           color[i] = colorData.nColors++;
                                                       } else {
                                                           colorData.nNodes[chosen_color]++;
                                                           color[i] = chosen_color;
                                                       }
                                                       
                                                       for (uint32_t j : neighborhood) mutex[j].unlock();
                                                   }
                                                   return colorData;
                                               },
                                               [](ColorData c1, ColorData c2) -> ColorData {
                                                   ColorData result;
                                                   result.nColors = std::max(c1.nColors, c2.nColors);
                                                   memset(result.nNodes, 0, sizeof(uint32_t) * result.nColors);
                                                   for (uint8_t i = 0; i < c1.nColors; ++i) result.nNodes[i] += c1.nNodes[i];
                                                   for (uint8_t i = 0; i < c2.nColors; ++i) result.nNodes[i] += c2.nNodes[i];
                                                   return result;
                                               });
    
    phases.resize(colorData.nColors);
    for (int i = 0; i < colorData.nColors; ++i) phases[i].reserve(colorData.nNodes[i]);
    
    for (uint32_t i = 0; i < size; ++i) phases[color[i]].push_back(i);
}
#else
void Hierarchy::generate_graph_coloring_deterministic(const AdjacentMatrix& adj, int size,
                                                      std::vector<std::vector<int>>& phases) {
    phases.clear();
    
    std::vector<uint32_t> perm(size);
    for (uint32_t i = 0; i < size; ++i) perm[i] = i;
    pcg32 rng;
    rng.shuffle(perm.begin(), perm.end());
    
    std::vector<int> color(size, -1);
    std::vector<uint8_t> possible_colors;
    std::vector<int> size_per_color;
    int ncolors = 0;
    
    for (uint32_t i = 0; i < size; ++i) {
        uint32_t ip = perm[i];
        
        std::fill(possible_colors.begin(), possible_colors.end(), 1);
        
        for (auto& link : adj[ip]) {
            int c = color[link.id];
            if (c >= 0) possible_colors[c] = 0;
        }
        
        int chosen_color = -1;
        for (uint32_t j = 0; j < possible_colors.size(); ++j) {
            if (possible_colors[j]) {
                chosen_color = j;
                break;
            }
        }
        
        if (chosen_color < 0) {
            chosen_color = ncolors++;
            possible_colors.resize(ncolors);
            size_per_color.push_back(0);
        }
        
        color[ip] = chosen_color;
        size_per_color[chosen_color]++;
    }
    phases.resize(ncolors);
    for (int i = 0; i < ncolors; ++i) phases[i].reserve(size_per_color[i]);
    for (uint32_t i = 0; i < size; ++i) phases[color[i]].push_back(i);
}
#endif

void Hierarchy::DownsampleGraph(const AdjacentMatrix adj, const MatrixXd& V, const MatrixXd& N,
                                const VectorXd& A, MatrixXd& V_p, MatrixXd& N_p, VectorXd& A_p,
                                MatrixXi& to_upper, VectorXi& to_lower, AdjacentMatrix& adj_p) {
    struct Entry {
        int i, j;
        double order;
        inline Entry() { i = j = -1; };
        inline Entry(int i, int j, double order) : i(i), j(j), order(order) {}
        inline bool operator<(const Entry& e) const { return order > e.order; }
        inline bool operator==(const Entry& e) const { return order == e.order; }
    };
    
    int nLinks = 0;
    for (auto& adj_i : adj) nLinks += adj_i.size();
    std::vector<Entry> entries(nLinks);
    std::vector<int> bases(adj.size());
    for (int i = 1; i < bases.size(); ++i) {
        bases[i] = bases[i - 1] + adj[i - 1].size();
    }
    
#ifdef WITH_OMP
#pragma omp parallel for
#endif
    for (int i = 0; i < V.cols(); ++i) {
        int base = bases[i];
        auto& ad = adj[i];
        auto entry_it = entries.begin() + base;
        for (auto it = ad.begin(); it != ad.end(); ++it, ++entry_it) {
            int k = it->id;
            double dp = N.col(i).dot(N.col(k));
            double ratio = A[i] > A[k] ? (A[i] / A[k]) : (A[k] / A[i]);
            *entry_it = Entry(i, k, dp * ratio);
        }
    }
    
    pss::parallel_stable_sort(entries.begin(), entries.end(), std::less<Entry>());
    
    std::vector<bool> mergeFlag(V.cols(), false);
    
    int nCollapsed = 0;
    for (int i = 0; i < nLinks; ++i) {
        const Entry& e = entries[i];
        if (mergeFlag[e.i] || mergeFlag[e.j]) continue;
        mergeFlag[e.i] = mergeFlag[e.j] = true;
        entries[nCollapsed++] = entries[i];
    }
    int vertexCount = V.cols() - nCollapsed;
    
    // Allocate memory for coarsened graph
    V_p.resize(3, vertexCount);
    N_p.resize(3, vertexCount);
    A_p.resize(vertexCount);
    to_upper.resize(2, vertexCount);
    to_lower.resize(V.cols());
    
#ifdef WITH_OMP
#pragma omp parallel for
#endif
    for (int i = 0; i < nCollapsed; ++i) {
        const Entry& e = entries[i];
        const double area1 = A[e.i], area2 = A[e.j], surfaceArea = area1 + area2;
        if (surfaceArea > RCPOVERFLOW)
            V_p.col(i) = (V.col(e.i) * area1 + V.col(e.j) * area2) / surfaceArea;
        else
            V_p.col(i) = (V.col(e.i) + V.col(e.j)) * 0.5f;
        Vector3d normal = N.col(e.i) * area1 + N.col(e.j) * area2;
        double norm = normal.norm();
        N_p.col(i) = norm > RCPOVERFLOW ? Vector3d(normal / norm) : Vector3d::UnitX();
        A_p[i] = surfaceArea;
        to_upper.col(i) << e.i, e.j;
        to_lower[e.i] = i;
        to_lower[e.j] = i;
    }
    
    int offset = nCollapsed;
    
    for (int i = 0; i < V.cols(); ++i) {
        if (!mergeFlag[i]) {
            int idx = offset++;
            V_p.col(idx) = V.col(i);
            N_p.col(idx) = N.col(i);
            A_p[idx] = A[i];
            to_upper.col(idx) << i, -1;
            to_lower[i] = idx;
        }
    }
    
    adj_p.resize(V_p.cols());
    std::vector<int> capacity(V_p.cols());
    std::vector<std::vector<Link>> scratches(V_p.cols());
#ifdef WITH_OMP
#pragma omp parallel for
#endif
    for (int i = 0; i < V_p.cols(); ++i) {
        int t = 0;
        for (int j = 0; j < 2; ++j) {
            int upper = to_upper(j, i);
            if (upper == -1) continue;
            t += adj[upper].size();
        }
        scratches[i].reserve(t);
        adj_p[i].reserve(t);
    }
#ifdef WITH_OMP
#pragma omp parallel for
#endif
    for (int i = 0; i < V_p.cols(); ++i) {
        auto& scratch = scratches[i];
        for (int j = 0; j < 2; ++j) {
            int upper = to_upper(j, i);
            if (upper == -1) continue;
            auto& ad = adj[upper];
            for (auto& link : ad) scratch.push_back(Link(to_lower[link.id], link.weight));
        }
        std::sort(scratch.begin(), scratch.end());
        int id = -1;
        auto& ad = adj_p[i];
        for (auto& link : scratch) {
            if (link.id != i) {
                if (id != link.id) {
                    ad.push_back(link);
                    id = link.id;
                } else {
                    ad.back().weight += link.weight;
                }
            }
        }
    }
}

void Hierarchy::SaveToFile(FILE* fp) {
    Save(fp, mScale);
    Save(fp, mF);
    Save(fp, mE2E);
    Save(fp, mAdj);
    Save(fp, mV);
    Save(fp, mN);
    Save(fp, mA);
    Save(fp, mToLower);
    Save(fp, mToUpper);
    Save(fp, mQ);
    Save(fp, mO);
    Save(fp, mS);
    Save(fp, mK);
    Save(fp, this->mPhases);
}

void Hierarchy::LoadFromFile(FILE* fp) {
    Read(fp, mScale);
    Read(fp, mF);
    Read(fp, mE2E);
    Read(fp, mAdj);
    Read(fp, mV);
    Read(fp, mN);
    Read(fp, mA);
    Read(fp, mToLower);
    Read(fp, mToUpper);
    Read(fp, mQ);
    Read(fp, mO);
    Read(fp, mS);
    Read(fp, mK);
    Read(fp, this->mPhases);
}

void Hierarchy::UpdateGraphValue(std::vector<Vector3i>& FQ, std::vector<Vector3i>& F2E,
                                 std::vector<Vector2i>& edge_diff) {
    FQ = std::move(mFQ[0]);
    F2E = std::move(mF2E[0]);
    edge_diff = std::move(mEdgeDiff[0]);
}

void Hierarchy::DownsampleEdgeGraph(std::vector<Vector3i>& FQ, std::vector<Vector3i>& F2E,
                                    std::vector<Vector2i>& edge_diff, int level) {
    std::vector<Vector2i> E2F(edge_diff.size(), Vector2i(-1, -1));
    for (int i = 0; i < F2E.size(); ++i) {
        for (int j = 0; j < 3; ++j) {
            int e = F2E[i][j];
            if (E2F[e][0] == -1)
                E2F[e][0] = i;
            else
                E2F[e][1] = i;
        }
    }
    int levels = (level == -1) ? 100 : level;
    mFQ.resize(levels);
    mF2E.resize(levels);
    mE2F.resize(levels);
    mEdgeDiff.resize(levels);
    mSing.resize(levels);
    mToUpperEdges.resize(levels - 1);
    mToUpperOrients.resize(levels - 1);
    for (int i = 0; i < FQ.size(); ++i) {
        Vector2i diff(0, 0);
        for (int j = 0; j < 3; ++j) {
            diff += rshift90(edge_diff[F2E[i][j]], FQ[i][j]);
        }
        if (diff != Vector2i::Zero()) {
            mSing[0].push_back(i);
        }
    }
    
    mFQ[0] = std::move(FQ);
    mF2E[0] = std::move(F2E);
    mE2F[0] = std::move(E2F);
    mEdgeDiff[0] = std::move(edge_diff);
    for (int l = 0; l < levels - 1; ++l) {
        auto& FQ = mFQ[l];
        auto& E2F = mE2F[l];
        auto& F2E = mF2E[l];
        auto& EdgeDiff = mEdgeDiff[l];
        auto& Sing = mSing[l];
        std::vector<int> fixed_faces(F2E.size(), 0);
        for (auto& s : Sing) {
            fixed_faces[s] = 1;
        }
        
        auto& toUpper = mToUpperEdges[l];
        auto& toUpperOrients = mToUpperOrients[l];
        toUpper.resize(E2F.size(), -1);
        toUpperOrients.resize(E2F.size(), 0);
        
        auto& nFQ = mFQ[l + 1];
        auto& nE2F = mE2F[l + 1];
        auto& nF2E = mF2E[l + 1];
        auto& nEdgeDiff = mEdgeDiff[l + 1];
        auto& nSing = mSing[l + 1];
        
        for (int i = 0; i < E2F.size(); ++i) {
            if (EdgeDiff[i] != Vector2i::Zero()) continue;
            if (fixed_faces[E2F[i][0]] || fixed_faces[E2F[i][1]]) {
                continue;
            }
            for (int j = 0; j < 2; ++j) {
                int f = E2F[i][j];
                for (int k = 0; k < 3; ++k) {
                    int neighbor_e = F2E[f][k];
                    for (int m = 0; m < 2; ++m) {
                        int neighbor_f = E2F[neighbor_e][m];
                        if (fixed_faces[neighbor_f] == 0) fixed_faces[neighbor_f] = 1;
                    }
                }
            }
            fixed_faces[E2F[i][0]] = 2;
            fixed_faces[E2F[i][1]] = 2;
            toUpper[i] = -2;
        }
        for (int i = 0; i < E2F.size(); ++i) {
            if (toUpper[i] == -2) continue;
            if (fixed_faces[E2F[i][0]] == 2 && fixed_faces[E2F[i][1]] == 2) {
                toUpper[i] = -3;
                continue;
            }
        }
        int numE = 0;
        for (int i = 0; i < toUpper.size(); ++i) {
            if (toUpper[i] == -1) {
                if (fixed_faces[E2F[i][0]] < 2 && fixed_faces[E2F[i][1]] < 2) {
                    nE2F.push_back(E2F[i]);
                    toUpperOrients[i] = 0;
                    toUpper[i] = numE++;
                    continue;
                }
                int f0 = fixed_faces[E2F[i][0]] < 2 ? E2F[i][0] : E2F[i][1];
                int e = i;
                int f = f0;
                std::vector<std::pair<int, int>> paths;
                paths.push_back(std::make_pair(i, 0));
                while (true) {
                    if (E2F[e][0] == f)
                        f = E2F[e][1];
                    else if (E2F[e][1] == f)
                        f = E2F[e][0];
                    if (fixed_faces[f] < 2) {
                        for (int j = 0; j < paths.size(); ++j) {
                            auto& p = paths[j];
                            toUpper[p.first] = numE;
                            int orient = p.second;
                            if (j > 0) orient = (orient + toUpperOrients[paths[j - 1].first]) % 4;
                            toUpperOrients[p.first] = orient;
                        }
                        nE2F.push_back(Vector2i(f0, f));
                        numE += 1;
                        break;
                    }
                    int ind0 = -1, ind1 = -1;
                    int e0 = e;
                    for (int j = 0; j < 3; ++j) {
                        if (F2E[f][j] == e) {
                            ind0 = j;
                            break;
                        }
                    }
                    for (int j = 0; j < 3; ++j) {
                        int e1 = F2E[f][j];
                        if (e1 != e && toUpper[e1] != -2) {
                            e = e1;
                            ind1 = j;
                            break;
                        }
                    }
                    
                    if (ind1 != -1) {
                        paths.push_back(std::make_pair(e, (FQ[f][ind1] - FQ[f][ind0] + 6) % 4));
                    } else {
                        if (EdgeDiff[e] != Vector2i::Zero()) {
                            printf("Unsatisfied !!!...\n");
                            printf("%d %d %d: %d %d\n", F2E[f][0], F2E[f][1], F2E[f][2], e0, e);
                            exit(0);
                        }
                        for (auto& p : paths) {
                            toUpper[p.first] = numE;
                            toUpperOrients[p.first] = 0;
                        }
                        numE += 1;
                        nE2F.push_back(Vector2i(f0, f0));
                        break;
                    }
                }
            }
        }
        nEdgeDiff.resize(numE);
        for (int i = 0; i < toUpper.size(); ++i) {
            if (toUpper[i] >= 0 && toUpperOrients[i] == 0) {
                nEdgeDiff[toUpper[i]] = EdgeDiff[i];
            }
        }
        std::vector<int> upperface(F2E.size(), -1);
        
        for (int i = 0; i < F2E.size(); ++i) {
            Vector3i eid;
            for (int j = 0; j < 3; ++j) {
                eid[j] = toUpper[F2E[i][j]];
            }
            if (eid[0] >= 0 && eid[1] >= 0 && eid[2] >= 0) {
                Vector3i eid_orient;
                for (int j = 0; j < 3; ++j) {
                    eid_orient[j] = (FQ[i][j] + 4 - toUpperOrients[F2E[i][j]]) % 4;
                }
                upperface[i] = nF2E.size();
                nF2E.push_back(eid);
                nFQ.push_back(eid_orient);
            }
        }
        for (int i = 0; i < nE2F.size(); ++i) {
            for (int j = 0; j < 2; ++j) {
                nE2F[i][j] = upperface[nE2F[i][j]];
            }
        }
        
        for (auto& s : Sing) {
            if (upperface[s] >= 0) nSing.push_back(upperface[s]);
        }
        mToUpperFaces.push_back(std::move(upperface));
        if (nEdgeDiff.size() == EdgeDiff.size()) {
            levels = l + 1;
            break;
        }
    }
    
    mFQ.resize(levels);
    mF2E.resize(levels);
    mE2F.resize(levels);
    mEdgeDiff.resize(levels);
    mSing.resize(levels);
    mToUpperEdges.resize(levels - 1);
    mToUpperOrients.resize(levels - 1);
}

void Hierarchy::FixFlip() {
    int l = mF2E.size() - 1;
    auto& F2E = mF2E[l];
    auto& E2F = mE2F[l];
    auto& FQ = mFQ[l];
    auto& EdgeDiff = mEdgeDiff[l];
    std::set<int> singular_edges;
    std::set<int> singular_faces;
    std::vector<int> E2E(F2E.size() * 3, -1);
    for (int i = 0; i < E2F.size(); ++i) {
        int v1 = E2F[i][0];
        int v2 = E2F[i][1];
        int t1 = 0;
        int t2 = 2;
        while (F2E[v1][t1] != i) t1 += 1;
        while (F2E[v2][t2] != i) t2 -= 1;
        t1 += v1 * 3;
        t2 += v2 * 3;
        E2E[t1] = t2;
        E2E[t2] = t1;
    }
    auto Area = [&](int f) {
        Vector2i diff1 = rshift90(EdgeDiff[F2E[f][0]], FQ[f][0]);
        Vector2i diff2 = rshift90(EdgeDiff[F2E[f][1]], FQ[f][1]);
        return diff1[0] * diff2[1] - diff1[1] * diff2[0];
    };
    std::vector<int> valences(F2E.size() * 3, -10000);
    int max_v = 0;
    int sings = 0;
    auto ComputeValence = [&](int i) {
        double angle = 0;
        int deid = i;
        do {
            int deid1 = deid / 3 * 3 + (deid + 2) % 3;
            auto diff1 = rshift90(EdgeDiff[F2E[deid/3][deid%3]], FQ[deid/3][deid%3]);
            auto diff2 = -rshift90(EdgeDiff[F2E[deid1/3][deid1%3]], FQ[deid1/3][deid1%3]);
            double a = atan2(diff1[0] * diff2[1] - diff1[1] * diff2[0], diff1[0] * diff2[0] + diff1[1] * diff2[1]);
            angle += a / 3.141592654 * 180;
            deid = E2E[deid1];
        } while (deid != i);
        return (int)((angle + 7210) / 90 - 80);
    };
    for (int i = 0; i < valences.size(); ++i) {
        if (valences[i] > -10000)
            continue;
        int deid = i;
        int v = ComputeValence(i);
        deid = i;
        do {
            valences[deid] = v;
            int deid1 = deid / 3 * 3 + (deid + 2) % 3;
            deid = E2E[deid1];
        } while (deid != i);
        max_v = std::max(max_v, v);
        if (v % 4 != 0) {
            sings += 1;
        }
    }
    auto CheckShrink = [&] (int deid, int allowed_edge) {
        if (deid == -1) {
            return false;
        }
        std::vector<int> corresponding_faces;
        std::vector<int> corresponding_edges;
        std::vector<Vector2i> corresponding_diff;
        Vector2i diff = EdgeDiff[F2E[deid / 3][deid % 3]];
        do {
            corresponding_diff.push_back(diff);
            corresponding_edges.push_back(deid);
            corresponding_faces.push_back(deid / 3);
            if (corresponding_faces.size() > F2E.size() * 10) {
                for (int i = 0; i < corresponding_edges.size(); ++i) {
                }
                exit(0);
            }
            // transform to the next face
            deid = E2E[deid];
            if (deid == -1) {
                return false;
            }
            // transform for the target incremental diff
            diff = -rshift90(diff, FQ[deid / 3][deid % 3]);
            deid = deid / 3 * 3 + (deid + 1) % 3;
            // transform to local
            diff = rshift90(diff, (4 - FQ[deid / 3][deid % 3]) % 4);
        } while (deid != corresponding_edges.front());
        // check diff
        if (diff != corresponding_diff.front()) {
            return false;
        }
        std::unordered_map<int, Vector2i> new_values;
        for (int i = 0; i < corresponding_diff.size(); ++i) {
            int deid = corresponding_edges[i];
            int eid = F2E[deid / 3][deid % 3];
            new_values[eid] = EdgeDiff[eid];
        }
        auto diffs = corresponding_diff;
        for (int i = 0; i < corresponding_diff.size(); ++i) {
            int deid = corresponding_edges[i];
            int eid = F2E[deid / 3][deid % 3];
            auto& res = new_values[eid];
            res -= corresponding_diff[i];
            int edge_thres = allowed_edge;
            if (abs(res[0]) > edge_thres ||
                abs(res[1]) > edge_thres) {
                return false;
            }
            if ((abs(res[0]) > 1 && abs(res[1]) != 0) ||
                (abs(res[1]) > 1 && abs(res[0]) != 0))
                return false;
        }
        int prev_area = 0, current_area = 0;
        for (int f = 0; f < corresponding_faces.size(); ++f) {
            int area = Area(corresponding_faces[f]);
            if (area < 0) prev_area += 1;
        }
        for (auto& p : new_values) {
            std::swap(EdgeDiff[p.first], p.second);
        }
        for (int f = 0; f < corresponding_faces.size(); ++f) {
            int area = Area(corresponding_faces[f]);
            if (area < 0) {
                current_area += 1;
            }
        }
        if (current_area < prev_area) {
            return true;
        }
        for (auto& p : new_values) {
            std::swap(EdgeDiff[p.first], p.second);
        }
        return false;
    };
    
    std::vector<int> eraseF(F2E.size(), 0);
    std::vector<int> AddValence(F2E.size(), 0);
    auto ShrinkQuad = [&](int deid) {
        int rdeid = E2E[deid];
        int d1 = E2E[deid / 3 * 3 + (deid + 1) % 3];
        int d2 = E2E[deid / 3 * 3 + (deid + 2) % 3];
        int r1 = E2E[rdeid / 3 * 3 + (rdeid + 1) % 3];
        int r2 = E2E[rdeid / 3 * 3 + (rdeid + 2) % 3];
        E2E[d1] = r2;
        E2E[r2] = d1;
        E2E[d2] = r1;
        E2E[r1] = d2;
        eraseF[deid / 3] = 1;
        eraseF[rdeid / 3] = 1;
        Vector2i diff_d1_origin = rshift90(EdgeDiff[F2E[d1/3][d1 % 3]], FQ[d1/3][d1%3]);
        Vector2i diff_r2 = EdgeDiff[F2E[r2/3][r2%3]];
        F2E[d1/3][d1%3] = F2E[r2/3][r2%3];
        int orient = 0;
        while (rshift90(diff_r2, orient) != diff_d1_origin)
            orient += 1;
        FQ[d1/3][d1%3] = orient;
        Vector2i diff_d2_origin = rshift90(EdgeDiff[F2E[d2/3][d2%3]], FQ[d2/3][d2%3]);
        Vector2i diff_r1 = EdgeDiff[F2E[r1/3][r1%3]];
        F2E[d2/3][d2%3] = F2E[r1/3][r1%3];
        orient = 0;
        while (rshift90(diff_r1, orient) != diff_d2_origin)
            orient += 1;
        FQ[d2/3][d2%3] = orient;
    };
    auto RemoveFlip = [&](int fid) {
        if (Area(fid) >= 0) {
            printf("Fail1\n");
            return false;
        }
        int deid = fid * 3;
        do {
            auto& diff = EdgeDiff[F2E[deid/3][deid%3]];
            if (abs(diff[0]) == 1 && abs(diff[1]) == 1)
                break;
            deid += 1;
        } while (true);
        // case one, the reverse is positive
        int rdeid = E2E[deid];
        if (Area(rdeid/3) > 0) {
            ShrinkQuad(deid);
            return true;
        } else {
            for (int d = 1; d < 3; ++d) {
                int rdeid = E2E[deid / 3 * 3 + (deid + d) % 3];
                if (Area(rdeid/3) > 0) {
                    int next_r = rdeid / 3 * 3;
                    do {
                        auto& diff = EdgeDiff[F2E[next_r/3][next_r%3]];
                        if (abs(diff[0]) == 1 && abs(diff[1]) == 1)
                            break;
                        next_r += 1;
                    } while (true);
                    if (Area(E2E[next_r]/3) > 0) {
                        ShrinkQuad(next_r);
                        ShrinkQuad(deid);
                        return true;
                    }
                }
            }
        }
        return false;
    };
    std::queue<int> flipped;
    for (int i = 0; i < F2E.size(); ++i) {
        int area = Area(i);
        if (area < 0) {
            flipped.push(i);
        }
    }
    
    bool update = false;
    int max_len = 1;
    while (!update && max_len < 2) {
        while (!flipped.empty()) {
            int f = flipped.front();
            if (Area(f) >= 0) {
                flipped.pop();
                continue;
            }
            for (int i = 0; i < 3; ++i) {
                if (CheckShrink(f * 3 + i, max_len) ||
                    CheckShrink(E2E[f * 3 + i], max_len)) {
                    update = true;
                    break;
                }
            }
            flipped.pop();
        }
        max_len += 1;
    }
    if (update) {
        Hierarchy flip_hierarchy;
        flip_hierarchy.DownsampleEdgeGraph(mFQ.back(), mF2E.back(), mEdgeDiff.back(), -1);
        flip_hierarchy.FixFlip();
        flip_hierarchy.UpdateGraphValue(mFQ.back(), mF2E.back(), mEdgeDiff.back());
    } else {
        /*
        update = false;
        std::queue<int> flipped, flipped_temp;
        for (int i = 0; i < F2E.size(); ++i) {
            int area = Area(i);
            if (area < 0) {
                flipped.push(i);
            }
        }
        while (flipped.size() > 0) {
            int flip_num = flipped.size();
            printf("Begin Flip: %d\n", flip_num);
            std::set<int> flipped_set;
            while (!flipped.empty()) {
                int f = flipped.front();
                flipped.pop();
                if (eraseF[f] || Area(f) >= 0 || RemoveFlip(f)) {
                    continue;
                }
                flipped_set.insert(f);
                flipped_temp.push(f);
            }
            std::swap(flipped, flipped_temp);
            printf("Final Flip: %d\n", flipped.size());
            if (flipped.size() == flip_num) {
                printf("Flip %d\n", flip_num);
                while (!flipped.empty()) {
                    int fid = flipped.front();
                    int deid = fid * 3;
                    do {
                        auto& diff = EdgeDiff[F2E[deid/3][deid%3]];
                        if (abs(diff[0]) == 1 && abs(diff[1]) == 1)
                            break;
                        deid += 1;
                    } while (true);
                    // case one, the reverse is positive
                    int rdeid = E2E[deid];
                    std::vector<int> fids;
                    fids.push_back(fid);
                    int fid1 = rdeid / 3;
                    fids.push_back(fid1);
                    for (int d = 1; d < 3; ++d) {
                        int rdeid = E2E[deid / 3 * 3 + (deid + d) % 3];
                        fids.push_back(rdeid/3);
                        if (Area(rdeid/3) > 0) {
                            int next_r = rdeid / 3 * 3;
                            do {
                                auto& diff = EdgeDiff[F2E[next_r/3][next_r%3]];
                                if (abs(diff[0]) == 1 && abs(diff[1]) == 1)
                                    break;
                                next_r += 1;
                            } while (true);
                            fids.push_back(E2E[next_r]/3);
                            if (Area(E2E[next_r]/3) > 0) {
                                ShrinkQuad(next_r);
                                ShrinkQuad(deid);
                            }
                        } else {
                            fids.push_back(-1);
                        }
                    }
                    for (auto& f : fids) {
                        printf("%d ", f);
                    }
                    printf("\n");
                    flipped.pop();
                }
                break;
            }
        }
        for (int i = 0; i < F2E.size(); ++i) {
            if (eraseF[i])
                continue;
            Vector2i diff(0, 0);
            for (int j = 0; j < 3; ++j) {
                diff += rshift90(EdgeDiff[F2E[i][j]], FQ[i][j]);
            }
            if (diff != Vector2i::Zero()) {
                printf("Non zero!\n");
            }
            if (Area(i) < 0) {
                printf("Flipp! %d %d %d\n", F2E[i][0], F2E[i][1], F2E[i][2]);
            }
            if (Area(i) == 0) {
                printf("Zero face!\n");
            }
        }
         */
    }
    PropagateEdge();
}

void Hierarchy::PropagateEdge()
{
    for (int level = mToUpperEdges.size(); level > 0; --level) {
        auto& EdgeDiff = mEdgeDiff[level];
        auto& nEdgeDiff = mEdgeDiff[level - 1];
        auto& FQ = mFQ[level];
        auto& nFQ = mFQ[level - 1];
        auto& F2E = mF2E[level - 1];
        auto& toUpper = mToUpperEdges[level - 1];
        auto& toUpperFace = mToUpperFaces[level - 1];
        auto& toUpperOrients = mToUpperOrients[level - 1];
        for (int i = 0; i < toUpper.size(); ++i) {
            if (toUpper[i] >= 0) {
                int orient = (4 - toUpperOrients[i]) % 4;
                nEdgeDiff[i] = rshift90(EdgeDiff[toUpper[i]], orient);
            } else {
                nEdgeDiff[i] = Vector2i(0, 0);
            }
        }
        for (int i = 0; i < toUpperFace.size(); ++i) {
            if (toUpperFace[i] == -1)
                continue;
            Vector3i eid_orient = FQ[toUpperFace[i]];
            for (int j = 0; j < 3; ++j) {
                nFQ[i][j] = (eid_orient[j] + toUpperOrients[F2E[i][j]]) % 4;
            }
        }
    }
}

#ifdef WITH_CUDA
#include <cuda_runtime.h>

void Hierarchy::CopyToDevice() {
    if (cudaAdj.empty()) {
        cudaAdj.resize(mAdj.size());
        cudaAdjOffset.resize(mAdj.size());
        for (int i = 0; i < mAdj.size(); ++i) {
            std::vector<int> offset(mAdj[i].size() + 1, 0);
            for (int j = 0; j < mAdj[i].size(); ++j) {
                offset[j + 1] = offset[j] + mAdj[i][j].size();
            }
            cudaMalloc(&cudaAdjOffset[i], sizeof(int) * (mAdj[i].size() + 1));
            cudaMemcpy(cudaAdjOffset[i], offset.data(), sizeof(int) * (mAdj[i].size() + 1),
                       cudaMemcpyHostToDevice);
            //            cudaAdjOffset[i] = (int*)malloc(sizeof(int) * (mAdj[i].size() + 1));
            //            memcpy(cudaAdjOffset[i], offset.data(), sizeof(int) * (mAdj[i].size() +
            //            1));
            
            cudaMalloc(&cudaAdj[i], sizeof(Link) * offset.back());
            //            cudaAdj[i] = (Link*)malloc(sizeof(Link) * offset.back());
            std::vector<Link> plainlink(offset.back());
            for (int j = 0; j < mAdj[i].size(); ++j) {
                memcpy(plainlink.data() + offset[j], mAdj[i][j].data(),
                       mAdj[i][j].size() * sizeof(Link));
            }
            cudaMemcpy(cudaAdj[i], plainlink.data(), plainlink.size() * sizeof(Link),
                       cudaMemcpyHostToDevice);
        }
    }
    
    if (cudaN.empty()) {
        cudaN.resize(mN.size());
        for (int i = 0; i < mN.size(); ++i) {
            cudaMalloc(&cudaN[i], sizeof(glm::dvec3) * mN[i].cols());
            //            cudaN[i] = (glm::dvec3*)malloc(sizeof(glm::dvec3) * mN[i].cols());
        }
    }
    for (int i = 0; i < mN.size(); ++i) {
        cudaMemcpy(cudaN[i], mN[i].data(), sizeof(glm::dvec3) * mN[i].cols(),
                   cudaMemcpyHostToDevice);
        //        memcpy(cudaN[i], mN[i].data(), sizeof(glm::dvec3) * mN[i].cols());
    }
    
    if (cudaV.empty()) {
        cudaV.resize(mV.size());
        for (int i = 0; i < mV.size(); ++i) {
            cudaMalloc(&cudaV[i], sizeof(glm::dvec3) * mV[i].cols());
            //            cudaV[i] = (glm::dvec3*)malloc(sizeof(glm::dvec3) * mV[i].cols());
        }
    }
    for (int i = 0; i < mV.size(); ++i) {
        cudaMemcpy(cudaV[i], mV[i].data(), sizeof(glm::dvec3) * mV[i].cols(),
                   cudaMemcpyHostToDevice);
        //        memcpy(cudaV[i], mV[i].data(), sizeof(glm::dvec3) * mV[i].cols());
    }
    
    if (cudaQ.empty()) {
        cudaQ.resize(mQ.size());
        for (int i = 0; i < mQ.size(); ++i) {
            cudaMalloc(&cudaQ[i], sizeof(glm::dvec3) * mQ[i].cols());
            //            cudaQ[i] = (glm::dvec3*)malloc(sizeof(glm::dvec3) * mQ[i].cols());
        }
    }
    for (int i = 0; i < mQ.size(); ++i) {
        cudaMemcpy(cudaQ[i], mQ[i].data(), sizeof(glm::dvec3) * mQ[i].cols(),
                   cudaMemcpyHostToDevice);
        //        memcpy(cudaQ[i], mQ[i].data(), sizeof(glm::dvec3) * mQ[i].cols());
    }
    if (cudaO.empty()) {
        cudaO.resize(mO.size());
        for (int i = 0; i < mO.size(); ++i) {
            cudaMalloc(&cudaO[i], sizeof(glm::dvec3) * mO[i].cols());
            //            cudaO[i] = (glm::dvec3*)malloc(sizeof(glm::dvec3) * mO[i].cols());
        }
    }
    for (int i = 0; i < mO.size(); ++i) {
        cudaMemcpy(cudaO[i], mO[i].data(), sizeof(glm::dvec3) * mO[i].cols(),
                   cudaMemcpyHostToDevice);
        //        memcpy(cudaO[i], mO[i].data(), sizeof(glm::dvec3) * mO[i].cols());
    }
    if (cudaPhases.empty()) {
        cudaPhases.resize(mPhases.size());
        for (int i = 0; i < mPhases.size(); ++i) {
            cudaPhases[i].resize(mPhases[i].size());
            for (int j = 0; j < mPhases[i].size(); ++j) {
                cudaMalloc(&cudaPhases[i][j], sizeof(int) * mPhases[i][j].size());
                //                cudaPhases[i][j] = (int*)malloc(sizeof(int) *
                //                mPhases[i][j].size());
            }
        }
    }
    for (int i = 0; i < mPhases.size(); ++i) {
        for (int j = 0; j < mPhases[i].size(); ++j) {
            cudaMemcpy(cudaPhases[i][j], mPhases[i][j].data(), sizeof(int) * mPhases[i][j].size(),
                       cudaMemcpyHostToDevice);
            //            memcpy(cudaPhases[i][j], mPhases[i][j].data(), sizeof(int) *
            //            mPhases[i][j].size());
        }
    }
    if (cudaToUpper.empty()) {
        cudaToUpper.resize(mToUpper.size());
        for (int i = 0; i < mToUpper.size(); ++i) {
            cudaMalloc(&cudaToUpper[i], mToUpper[i].cols() * sizeof(glm::ivec2));
            //            cudaToUpper[i] = (glm::ivec2*)malloc(mToUpper[i].cols() *
            //            sizeof(glm::ivec2));
        }
    }
    for (int i = 0; i < mToUpper.size(); ++i) {
        cudaMemcpy(cudaToUpper[i], mToUpper[i].data(), sizeof(glm::ivec2) * mToUpper[i].cols(),
                   cudaMemcpyHostToDevice);
        //        memcpy(cudaToUpper[i], mToUpper[i].data(), sizeof(glm::ivec2) *
        //        mToUpper[i].cols());
    }
    cudaDeviceSynchronize();
}

void Hierarchy::CopyToHost() {}

#endif

