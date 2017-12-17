#include "dedge.h"

#include "CompareKey.h"

#include <fstream>

void compute_direct_graph(MatrixXf& V, MatrixXi& F, VectorXi& V2E,
	VectorXi& E2E, VectorXi& boundary, VectorXi& nonManifold)
{
	printf("stage 1...\n");
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
	printf("stage 2...\n");

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
	printf("stage 3...\n");

	for (int i = 0; i < V2E.rows(); ++i) {
		int idx = V2E[i];
		int start = idx;
		int v2e = INT_MAX;
		do {
			if (i == 24)
				printf("%d\n", idx);
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

