#include "config.hpp"
#include "optimizer.hpp"
#include "field-math.hpp"
#include "flow.hpp"
#include <fstream>
#include <Eigen/Sparse>
#ifdef WITH_CUDA
#include <cuda_runtime.h>
#endif
#include "parametrizer.hpp"
Optimizer::Optimizer()
{}

void Optimizer::optimize_orientations(Hierarchy &mRes)
{
#ifdef WITH_CUDA
	optimize_orientations_cuda(mRes);
	printf("%s\n", cudaGetErrorString(cudaDeviceSynchronize()));
	cudaMemcpy(mRes.mQ[0].data(), mRes.cudaQ[0], sizeof(glm::dvec3) * mRes.mQ[0].cols(), cudaMemcpyDeviceToHost);

#else
	
	int levelIterations = 6;
	for (int level = mRes.mN.size() - 1; level >= 0; --level) {
		AdjacentMatrix &adj = mRes.mAdj[level];
		const MatrixXd &N = mRes.mN[level];
		MatrixXd &Q = mRes.mQ[level];
		auto& phases = mRes.mPhases[level];
		for (int iter = 0; iter < levelIterations; ++iter) {
			for (int phase = 0; phase < phases.size(); ++phase) {
				auto& p = phases[phase];
#ifdef WITH_OMP
#pragma omp parallel for
#endif
				for (int pi = 0; pi < p.size(); ++pi) {
					int i = p[pi];
					const Vector3d n_i = N.col(i);
					double weight_sum = 0.0f;
					Vector3d sum = Q.col(i);
					for (auto& link : adj[i]) {
						const int j = link.id;
						const double weight = link.weight;
						if (weight == 0)
							continue;
						const Vector3d n_j = N.col(j);
						Vector3d q_j = Q.col(j);
						std::pair<Vector3d, Vector3d> value = compat_orientation_extrinsic_4(sum, n_i, q_j, n_j);
						sum = value.first * weight_sum + value.second * weight;
						sum -= n_i*n_i.dot(sum);
						weight_sum += weight;
						double norm = sum.norm();
						if (norm > RCPOVERFLOW)
							sum /= norm;
					}
					if (weight_sum > 0) {
						Q.col(i) = sum;
					}
				}
			}
		}
		if (level > 0) {
			const MatrixXd &srcField = mRes.mQ[level];
			const MatrixXi &toUpper = mRes.mToUpper[level - 1];
			MatrixXd &destField = mRes.mQ[level - 1];
			const MatrixXd &N = mRes.mN[level - 1];
#ifdef WITH_OMP
#pragma omp parallel for
#endif
			for (int i = 0; i < srcField.cols(); ++i) {
				for (int k = 0; k < 2; ++k) {
					int dest = toUpper(k, i);
					if (dest == -1)
						continue;
					Vector3d q = srcField.col(i), n = N.col(dest);
					destField.col(dest) = q - n * n.dot(q);
				}
			}
		}
	}

	for (int l = 0; l< mRes.mN.size() - 1; ++l)  {
		const MatrixXd &N = mRes.mN[l];
		const MatrixXd &N_next = mRes.mN[l + 1];
		const MatrixXd &Q = mRes.mQ[l];
		MatrixXd &Q_next = mRes.mQ[l + 1];
		auto& toUpper = mRes.mToUpper[l];
#ifdef WITH_OMP
#pragma omp parallel for
#endif
		for (int i = 0; i < toUpper.cols(); ++i) {
			Vector2i upper = toUpper.col(i);
			Vector3d q0 = Q.col(upper[0]);
			Vector3d n0 = N.col(upper[0]);
			Vector3d q;

			if (upper[1] != -1) {
				Vector3d q1 = Q.col(upper[1]);
				Vector3d n1 = N.col(upper[1]);
				auto result = compat_orientation_extrinsic_4(q0, n0, q1, n1);
				q = result.first + result.second;
			}
			else {
				q = q0;
			}
			Vector3d n = N_next.col(i);
			q -= n.dot(q) * n;
			if (q.squaredNorm() > RCPOVERFLOW)
				q.normalize();

			Q_next.col(i) = q;
		}
	}

#endif
}

void Optimizer::optimize_scale(Hierarchy &mRes)
{
	const MatrixXd &N = mRes.mN[0];
	MatrixXd &Q = mRes.mQ[0];
	MatrixXd &V = mRes.mV[0];
	MatrixXd &S = mRes.mS[0];
	MatrixXd &K = mRes.mK[0];
	MatrixXi &F = mRes.mF;
	std::vector<Eigen::Triplet<double> > lhsTriplets;

	lhsTriplets.reserve(F.cols() * 6);
	std::vector<std::map<int, double> > entries(V.cols() * 2);
	for (int i = 0; i < F.cols(); ++i) {
		for (int j = 0; j < 3; ++j) {
			int v1 = F(j, i);
			int v2 = F((j + 1) % 3, i);
			Vector3d diff = V.col(v2) - V.col(v1);
			Vector3d q_1 = Q.col(v1);
			Vector3d q_2 = Q.col(v2);
			Vector3d n_1 = N.col(v1);
			Vector3d n_2 = N.col(v2);
			Vector3d q_1_y = n_1.cross(q_1);
			auto index = compat_orientation_extrinsic_index_4(q_1, n_1, q_2, n_2);
			int v1_x = v1 * 2, v1_y = v1 * 2 + 1, v2_x = v2 * 2, v2_y = v2 * 2 + 1;

			double dx = diff.dot(q_1);
			double dy = diff.dot(q_1_y);

			double kx_g = K(0, v1);
			double ky_g = K(1, v1);

			if (index.first % 2 != index.second % 2) {
				std::swap(v2_x, v2_y);
			}
			double scale_x = log(fmin(fmax(1 + kx_g * dy, 0.1), 10));
			double scale_y = log(fmin(fmax(1 + ky_g * dx, 0.1), 10));

			auto it = entries[v1_x].find(v2_x);
			if (it == entries[v1_x].end()) {
				entries[v1_x][v2_x] = -scale_x;
				entries[v2_x][v1_x] = scale_x;
				entries[v1_y][v2_y] = -scale_y;
				entries[v2_y][v1_y] = scale_y;
			}
			else {
				it->second -= scale_x;
				entries[v2_x][v1_x] += scale_x;
				entries[v1_y][v2_y] -= scale_y;
				entries[v2_y][v1_y] += scale_y;
			}
		}
	}

	for (int i = 0; i < entries.size(); ++i) {
		for (auto& j : entries[i]) {
			if (abs(j.second + entries[j.first][i]) > 1e-6) {
				printf("not true... %f %f\n", j.second, entries[j.first][i]);
                exit(0);
			}
		}
	}

	Eigen::SparseMatrix<double> A(V.cols() * 2, V.cols() * 2);
	VectorXd rhs(V.cols() * 2);
	rhs.setZero();
	for (int i = 0; i < entries.size(); ++i) {
		lhsTriplets.push_back(Eigen::Triplet<double>(i, i, entries[i].size()));
		for (auto& rec : entries[i]) {
			rhs(i) += rec.second;
			lhsTriplets.push_back(Eigen::Triplet<double>(i, rec.first, -1));
		}
	}
	A.setFromTriplets(lhsTriplets.begin(), lhsTriplets.end());
	Eigen::SparseLU< Eigen::SparseMatrix<double> > solver;
	solver.analyzePattern(A);

	solver.factorize(A);

	VectorXd result = solver.solve(rhs);

	for (int i = 0; i < V.cols(); ++i) {
		S(0, i) = exp(result(i * 2));
		S(1, i) = exp(result(i * 2 + 1));
	}

	for (int l = 0; l< mRes.mS.size() - 1; ++l)  {
		const MatrixXd &S = mRes.mS[l];
		MatrixXd &S_next = mRes.mS[l + 1];
		auto& toUpper = mRes.mToUpper[l];
		for (int i = 0; i < toUpper.cols(); ++i) {
			Vector2i upper = toUpper.col(i);
			Vector2d q0 = S.col(upper[0]);

			if (upper[1] != -1) {
				q0 = (q0 + S.col(upper[1])) * 0.5;
			}
			S_next.col(i) = q0;
		}
	}
}

void Optimizer::optimize_positions(Hierarchy &mRes, int with_scale)
{
	int levelIterations = 6;
#ifdef WITH_CUDA
	optimize_positions_cuda(mRes);
	cudaMemcpy(mRes.mO[0].data(), mRes.cudaO[0], sizeof(glm::dvec3) * mRes.mO[0].cols(), cudaMemcpyDeviceToHost);
#else
	for (int level = mRes.mAdj.size() - 1; level >= 0; --level) {
		for (int iter = 0; iter < levelIterations; ++iter) {
			AdjacentMatrix &adj = mRes.mAdj[level];
			const MatrixXd &N = mRes.mN[level], &Q = mRes.mQ[level], &V = mRes.mV[level];
			MatrixXd &O = mRes.mO[level];
			auto& phases = mRes.mPhases[level];
			for (int phase = 0; phase < phases.size(); ++phase) {
				auto& p = phases[phase];
#ifdef WITH_OMP
#pragma omp parallel for
#endif
				for (int pi = 0; pi < p.size(); ++pi) {
					int i = p[pi];
					double scale_x = mRes.mScale;
					double scale_y = mRes.mScale;
//					if (with_scale) {
//						scale_x *= S(0, i);
//						scale_y *= S(1, i);
//					}
					double inv_scale_x = 1.0f / scale_x;
					double inv_scale_y = 1.0f / scale_y;
					const Vector3d n_i = N.col(i), v_i = V.col(i);
					Vector3d q_i = Q.col(i);

					Vector3d sum = O.col(i);
					double weight_sum = 0.0f;

					q_i.normalize();
					for (auto& link : adj[i]) {
						const int j = link.id;
						const double weight = link.weight;
						if (weight == 0)
							continue;
						double scale_x_1 = mRes.mScale;
						double scale_y_1 = mRes.mScale;
						if (with_scale) {
//							scale_x_1 *= S(0, j);
//							scale_y_1 *= S(1, j);
						}
						double inv_scale_x_1 = 1.0f / scale_x_1;
						double inv_scale_y_1 = 1.0f / scale_y_1;

						const Vector3d n_j = N.col(j), v_j = V.col(j);
						Vector3d q_j = Q.col(j), o_j = O.col(j);

						q_j.normalize();

						std::pair<Vector3d, Vector3d> value = compat_position_extrinsic_4(
							v_i, n_i, q_i, sum, v_j, n_j, q_j, o_j,
							scale_x, scale_y, inv_scale_x, inv_scale_y,
							scale_x_1, scale_y_1, inv_scale_x_1, inv_scale_y_1);

						sum = value.first*weight_sum + value.second*weight;
						weight_sum += weight;
						if (weight_sum > RCPOVERFLOW)
							sum /= weight_sum;
						sum -= n_i.dot(sum - v_i)*n_i;
					}

					if (weight_sum > 0) {
						O.col(i) = position_round_4(sum, q_i, n_i, v_i, scale_x, scale_y, inv_scale_x, inv_scale_y);
					}
				}
			}
		}
		if (level > 0) {
			const MatrixXd &srcField = mRes.mO[level];
			const MatrixXi &toUpper = mRes.mToUpper[level - 1];
			MatrixXd &destField = mRes.mO[level - 1];
			const MatrixXd &N = mRes.mN[level - 1];
			const MatrixXd &V = mRes.mV[level - 1];
#ifdef WITH_OMP
#pragma omp parallel for
#endif
			for (int i = 0; i < srcField.cols(); ++i) {
				for (int k = 0; k < 2; ++k) {
					int dest = toUpper(k, i);
					if (dest == -1)
						continue;
					Vector3d o = srcField.col(i), n = N.col(dest), v = V.col(dest);
					o -= n * n.dot(o - v);
					destField.col(dest) = o;
				}
			}
		}
	}
#endif
}

void Optimizer::optimize_integer_constraints(Hierarchy &mRes, std::map<int, int>& singularities)
{
    std::vector<std::set<int> > singular_edges(mRes.mF2E.size());
    for (auto& f : singularities) {
        for (int j = 0; j < 3; ++j)
            singular_edges[0].insert(mRes.mF2E[0][f.first][j]);
    }
    for (int level = 0; level < mRes.mToUpperEdges.size(); ++level) {
        auto& toUpper = mRes.mToUpperEdges[level];
        auto& SingEdges = singular_edges[level];
        auto& nSingEdges = singular_edges[level + 1];
        for (auto& e : SingEdges) {
            if (toUpper[e] >= 0)
                nSingEdges.insert(toUpper[e]);
        }
    }
    
    bool FullFlow = false;
    for (int level = mRes.mToUpperEdges.size(); level >= 0; --level) {
        printf("flow %d...\n", level);
        auto& EdgeDiff = mRes.mEdgeDiff[level];
        auto& FQ = mRes.mFQ[level];
        auto& F2E = mRes.mF2E[level];
        auto& E2F = mRes.mE2F[level];
        auto& SingEdges = singular_edges[level];
        if (!FullFlow) {
            std::vector<Vector4i> edge_to_constraints(E2F.size() * 2, Vector4i(-1, 0, -1, 0));
            std::vector<int> initial(F2E.size() * 2, 0);
            for (int i = 0; i < F2E.size(); ++i) {
                for (int j = 0; j < 3; ++j) {
                    int e = F2E[i][j];
                    Vector2i index = rshift90(Vector2i(e * 2 + 1, e * 2 + 2),
                                          FQ[i][j]);
                    for (int k = 0; k < 2; ++k) {
                        int l = abs(index[k]);
                        int s = index[k] / l;
                        int ind = l - 1;
                        int equationID = i * 2 + k;
                        if (edge_to_constraints[ind][0] == -1) {
                            edge_to_constraints[ind][0] = equationID;
                            edge_to_constraints[ind][1] = s;
                        } else {
                            edge_to_constraints[ind][2] = equationID;
                            edge_to_constraints[ind][3] = s;
                        }
                        initial[equationID] += s * EdgeDiff[ind/2][ind%2];
                    }
                }
            }
            std::vector<std::pair<Vector2i, int> > arcs;
            std::vector<int> arc_ids;
            std::vector<int> singularity_edge;
            for (int i = 0; i < edge_to_constraints.size(); ++i) {
                if (edge_to_constraints[i][1] == -edge_to_constraints[i][3]) {
                    int v1 = edge_to_constraints[i][0];
                    int v2 = edge_to_constraints[i][2];
                    if (edge_to_constraints[i][1] < 0)
                        std::swap(v1, v2);
                    int current_v = EdgeDiff[i / 2][i % 2];
                    arcs.push_back(std::make_pair(Vector2i(v1, v2), current_v));
                    singularity_edge.push_back(SingEdges.count(i / 2));
                    arc_ids.push_back(i);
                }
            }
            int supply = 0;
            int demand = 0;
            for (int i = 0; i < initial.size(); ++i) {
                int init_val = initial[i];
                if (init_val > 0) {
                    arcs.push_back(std::make_pair(Vector2i(-1, i), initial[i]));
                    supply += init_val;
                }
                else if (init_val < 0) {
                    demand -= init_val;
                    arcs.push_back(std::make_pair(Vector2i(i, initial.size()), -init_val));
                }
            }
            
            MaxFlowHelper* flow = 0;
            if (supply < 5)
                flow = new BoykovMaxFlowHelper;
            else
                flow = new BoykovMaxFlowHelper;

            //    flow.resize(constraints_index.size() + 2, edge_diff.size() * 2);
            flow->resize(initial.size() + 2, arc_ids.size());
            std::set<int> ids;
            for (int i = 0; i < arcs.size(); ++i) {
                int v1 = arcs[i].first[0] + 1;
                int v2 = arcs[i].first[1] + 1;
                int c = arcs[i].second;
                if (v1 == 0 || v2 == initial.size() + 1) {
                    flow->AddEdge(v1, v2, c, 0, -1);
                }
                else {
                    flow->AddEdge(v1, v2, std::max(0, c + 2), std::max(0, -c + 2), arc_ids[i]);
                }
            }
            
            printf("compute....\n");
            int flow_count = flow->compute();
            printf("finish...\n");
            //flow.compute(edge_diff, face_edgeIds, face_edgeOrients, true);
            //    flow_count += flow.compute(edge_diff, face_edgeIds, face_edgeOrients, false);
            flow->Apply(EdgeDiff);
            delete flow;
            printf("%d %d %d\n", flow_count, supply, demand);
            if (flow_count == supply) {
                FullFlow = true;
            }
        }
        if (level != 0) {
            auto& nEdgeDiff = mRes.mEdgeDiff[level - 1];
            auto& toUpper = mRes.mToUpperEdges[level - 1];
            auto& toUpperOrients = mRes.mToUpperOrients[level - 1];
            for (int i = 0; i < toUpper.size(); ++i) {
                if (toUpper[i] >= 0) {
                    int orient = (4 - toUpperOrients[i]) % 4;
                    nEdgeDiff[i] = rshift90(EdgeDiff[toUpper[i]], orient);
                }
            }
        }
        if (FullFlow) {
            auto& EdgeDiff = mRes.mEdgeDiff[level];
            auto& F2E = mRes.mF2E[level];
            auto& FQ = mRes.mFQ[level];
            for (int i = 0; i < FQ.size(); ++i) {
                Vector2i d(0, 0);
                for (int j = 0; j < 3; ++j) {
                    d += rshift90(EdgeDiff[F2E[i][j]], FQ[i][j]);
                }
                if (d != Vector2i::Zero()) {
                    printf("wrong..\n");
                    exit(0);
                }
            }
        }
    }
    printf("finish...\n");
}


#ifdef WITH_CUDA

void Optimizer::optimize_orientations_cuda(Hierarchy &mRes)
{
	int levelIterations = 6;
	for (int level = mRes.mN.size() - 1; level >= 0; --level) {
		Link* adj = mRes.cudaAdj[level];
		int* adjOffset = mRes.cudaAdjOffset[level];
		glm::dvec3* N = mRes.cudaN[level];
		glm::dvec3* Q = mRes.cudaQ[level];
		auto& phases = mRes.cudaPhases[level];
		for (int iter = 0; iter < levelIterations; ++iter) {
			for (int phase = 0; phase < phases.size(); ++phase) {
				int* p = phases[phase];
				UpdateOrientation(p, mRes.mPhases[level][phase].size(), N, Q, adj, adjOffset, mRes.mAdj[level][phase].size());
			}
		}
		if (level > 0) {
			glm::dvec3* srcField = mRes.cudaQ[level];
			glm::ivec2* toUpper = mRes.cudaToUpper[level - 1];
			glm::dvec3* destField = mRes.cudaQ[level - 1];
			glm::dvec3* N = mRes.cudaN[level - 1];
			PropagateOrientationUpper(srcField, mRes.mQ[level].cols(), toUpper, N, destField);
		}
	}

	for (int l = 0; l< mRes.mN.size() - 1; ++l)  {
		glm::dvec3* N = mRes.cudaN[l];
		glm::dvec3* N_next = mRes.cudaN[l + 1];
		glm::dvec3* Q = mRes.cudaQ[l];
		glm::dvec3* Q_next = mRes.cudaQ[l + 1];
		glm::ivec2* toUpper = mRes.cudaToUpper[l];

		PropagateOrientationLower(toUpper, Q, N, Q_next, N_next, mRes.mToUpper[l].cols());
	}
}

void Optimizer::optimize_positions_cuda(Hierarchy &mRes)
{
	int levelIterations = 6;
	for (int level = mRes.mAdj.size() - 1; level >= 0; --level) {
		Link* adj = mRes.cudaAdj[level];
		int* adjOffset = mRes.cudaAdjOffset[level];
		glm::dvec3* N = mRes.cudaN[level];
		glm::dvec3* Q = mRes.cudaQ[level];
		glm::dvec3* V = mRes.cudaV[level];
		glm::dvec3* O = mRes.cudaO[level];
		std::vector<int*> phases = mRes.cudaPhases[level];
		for (int iter = 0; iter < levelIterations; ++iter) {
			for (int phase = 0; phase < phases.size(); ++phase) {
				int* p = phases[phase];
				UpdatePosition(p, mRes.mPhases[level][phase].size(), N, Q, adj, adjOffset, mRes.mAdj[level][phase].size(),
					V, O, mRes.mScale);
			}
		}
		if (level > 0) {
			glm::dvec3* srcField = mRes.cudaO[level];
			glm::ivec2* toUpper = mRes.cudaToUpper[level - 1];
			glm::dvec3* destField = mRes.cudaO[level - 1];
			glm::dvec3* N = mRes.cudaN[level - 1];
			glm::dvec3* V = mRes.cudaV[level - 1];
			PropagatePositionUpper(srcField, mRes.mO[level].cols(), toUpper, N, V, destField);
		}
	}
}

#endif
