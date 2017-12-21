#include "dedge.h"

#include "CompareKey.h"
#include <vector>
#include <fstream>


inline int dedge_prev(int e, int deg) { return (e % deg == 0u) ? e + (deg - 1) : e - 1; }

void compute_direct_graph(MatrixXd& V, MatrixXi& F, VectorXi& V2E,
	VectorXi& E2E, VectorXi& boundary, VectorXi& nonManifold)
{
	V2E.resize(V.cols());
	V2E.setConstant(-1);

	unsigned int deg = F.rows();
	std::vector<std::pair<unsigned int, unsigned int>> tmp(F.size());

	for (int f = 0; f < F.cols(); ++f) {
		for (unsigned int i = 0; i < deg; ++i) {
			unsigned int idx_cur = F(i, f),
				idx_next = F((i + 1) % deg, f),
				edge_id = deg * f + i;
			if (idx_cur >= V.cols() || idx_next >= V.cols())
				throw std::runtime_error("Mesh data contains an out-of-bounds vertex reference!");
			if (idx_cur == idx_next)
				continue;

			tmp[edge_id] = std::make_pair(idx_next, -1);
			if (V2E[idx_cur] == -1)
				V2E[idx_cur] = edge_id;
			else {
				unsigned int idx = V2E[idx_cur];
				while (tmp[idx].second != -1) {
					idx = tmp[idx].second;
				}
				tmp[idx].second = edge_id;
			}
		}
	}

	nonManifold.resize(V.cols());
	nonManifold.setConstant(false);

	E2E.resize(F.cols() * deg);
	E2E.setConstant(-1);
	for (int f = 0; f < F.cols(); ++f) {
		for (unsigned int i = 0; i < deg; ++i) {
			unsigned int idx_cur = F(i, f),
				idx_next = F((i + 1) % deg, f),
				edge_id_cur = deg * f + i;
			if (idx_cur == idx_next)
				continue;

			unsigned int it = V2E[idx_next], edge_id_opp = -1;
			while (it != -1) {
				if (tmp[it].first == idx_cur) {
					if (edge_id_opp == -1) {
						edge_id_opp = it;
					}
					else {
						nonManifold[idx_cur] = true;
						nonManifold[idx_next] = true;
						edge_id_opp = -1;
						break;
					}
				}
				it = tmp[it].second;
			}

			if (edge_id_opp != -1 && edge_id_cur < edge_id_opp) {
				E2E[edge_id_cur] = edge_id_opp;
				E2E[edge_id_opp] = edge_id_cur;
			}
		}
	}

	unsigned int nonManifoldCounter(0), boundaryCounter(0), isolatedCounter(0);

	boundary.resize(V.cols());
	boundary.setConstant(false);

	/* Detect boundary regions of the mesh and adjust vertex->edge pointers*/
	for (int i = 0; i < V.cols(); ++i) {
		unsigned int edge = V2E[i];
		if (edge == -1) {
			isolatedCounter++;
			continue;
		}
		if (nonManifold[i]) {
			nonManifoldCounter++;
			V2E[i] = -1;
			continue;
		}

		/* Walk backwards to the first boundary edge (if any) */
		unsigned int start = edge, v2e = -1;
		do {
			v2e = std::min(v2e, edge);
			unsigned int prevEdge = E2E[dedge_prev(edge, deg)];
			if (prevEdge == -1) {
				/* Reached boundary -- update the vertex->edge link */
				v2e = edge;
				boundary[i] = true;
				boundaryCounter++;
				break;
			}
			edge = prevEdge;
		} while (edge != start);
		V2E[i] = v2e;
	}
}

/*

void compute_direct_graph(MatrixXd& V, MatrixXi& F, VectorXi& V2E,
VectorXi& E2E, VectorXi& boundary, VectorXi& nonManifold)
{
	E2E.resize(F.cols() * 3);
	std::map<Key2i, std::pair<int, int> > edge_indices;
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < F.cols(); ++j) {
			int id1 = F(i, j);
			int id2 = F((i + 1) % 3, j);
			auto edge_id = Key2i(id2, id1);
			if (edge_indices.count(edge_id) == 0) {
				edge_id = Key2i(id1, id2);
				edge_indices[edge_id] = std::make_pair(j * 3 + i, -1);
				E2E[j * 3 + i] = -1;
			}
			else {
				auto& p = edge_indices[edge_id];
				if (p.second == -1) {
					p.second = j * 3 + i;
					E2E[p.first] = p.second;
					E2E[p.second] = p.first;
				}
				else {
					E2E[p.first] = -2;
					E2E[p.second] = -2;
					E2E[j * 3 + i] = -2;
				}
			}
		}
	}

	V2E.resize(V.cols());
	V2E.setConstant(-1);
	boundary.resize(V.cols());
	boundary.setZero();
	nonManifold.resize(V.cols());
	nonManifold.setZero();
	for (int i = 0; i < E2E.rows(); ++i) {
		int row = i % 3;
		int col = i / 3;
		int id1 = F(row, col);
		int id2 = F((row + 1) % 3, col);
		if (E2E[i] == -2) {
			nonManifold[id1] = 1;
			nonManifold[id2] = 1;
		}
		else {
			if (V2E[id1] == -1)
				V2E[id1] = i;
			if (E2E[i] == -1) {
				boundary[id1] = 1;
				boundary[id2] = 1;
			}
		}
	}

	for (int i = 0; i < V2E.rows(); ++i) {
		int idx = V2E[i];
		int start = idx;
		int v2e = INT_MAX;
		do {
			if (idx < v2e)
				v2e = idx;
			int f_ind = idx % 3;
			int new_idx = E2E[idx - f_ind + (f_ind + 2) % 3];
			if (new_idx == -1) {
				v2e = idx;
				break;
			}
			idx = new_idx;
		} while (start != idx);
		V2E[i] = v2e;
	}
}

*/