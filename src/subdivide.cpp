#include "subdivide.hpp"
#include "dedge.hpp"
#include <queue>
#include <fstream>

void subdivide(MatrixXi &F, MatrixXd &V, VectorXi &V2E, VectorXi &E2E,
	VectorXi &boundary, VectorXi &nonmanifold, double maxLength) {
	typedef std::pair<double, int> Edge;

	std::priority_queue<Edge> queue;

	maxLength *= maxLength;

	for (int i = 0; i < E2E.size(); ++i) {
		int v0 = F(i % 3, i / 3), v1 = F((i + 1) % 3, i / 3);
		if (nonmanifold[v0] || nonmanifold[v1])
			continue;
		double length = (V.col(v0) - V.col(v1)).squaredNorm();
		if (length > maxLength) {
			int other = E2E[i];
			if (other == -1 || other > i)
				queue.push(Edge(length, i));
		}
	}

	int nV = V.cols(), nF = F.cols(), nSplit = 0;
	/*
	/   v0  \
	v1p 1 | 0 v0p
	\   v1  /

	/   v0  \
	/  1 | 0  \
	v1p - vn - v0p
	\  2 | 3  /
	\   v1  /

	f0: vn, v0p, v0
	f1: vn, v0, v1p
	f2: vn, v1p, v1
	f3: vn, v1, v0p
	*/
	int counter = 0;
	while (!queue.empty()) {
		counter += 1;
		Edge edge = queue.top();
		queue.pop();
		int e0 = edge.second, e1 = E2E[e0];
		bool is_boundary = e1 == -1;
		int f0 = e0 / 3, f1 = is_boundary ? -1 : (e1 / 3);
		int v0 = F(e0 % 3, f0), v0p = F((e0 + 2) % 3, f0), v1 = F((e0 + 1) % 3, f0);
		if ((V.col(v0) - V.col(v1)).squaredNorm() != edge.first) {
			continue;
		}
		int v1p = is_boundary ? -1 : F((e1 + 2) % 3, f1);
		int vn = nV++;
		nSplit++;
		/* Update V */
		if (nV > V.cols()) {
			V.conservativeResize(V.rows(), V.cols() * 2);
			V2E.conservativeResize(V.cols());
			boundary.conservativeResize(V.cols());
			nonmanifold.conservativeResize(V.cols());
		}

		/* Update V */
		V.col(vn) = (V.col(v0) + V.col(v1)) * 0.5f;
		nonmanifold[vn] = false;
		boundary[vn] = is_boundary;

		/* Update F and E2E */
		int f2 = is_boundary ? -1 : (nF++);
		int f3 = nF++;
		int f20, f30;
		if (nF > F.cols()) {
			F.conservativeResize(F.rows(), std::max(nF, (int)F.cols() * 2));
			E2E.conservativeResize(F.cols() * 3);
		}

		/* Update F */
		F.col(f0) << vn, v0p, v0;
		if (!is_boundary) {
			F.col(f1) << vn, v0, v1p;
			F.col(f2) << vn, v1p, v1;
		}
		F.col(f3) << vn, v1, v0p;

		/* Update E2E */
		const int e0p = E2E[dedge_prev_3(e0)],
			e0n = E2E[dedge_next_3(e0)];

#define sE2E(a, b) E2E[a] = b; if (b != -1) E2E[b] = a;
		sE2E(3 * f0 + 0, 3 * f3 + 2);
		sE2E(3 * f0 + 1, e0p);
		sE2E(3 * f3 + 1, e0n);
		if (is_boundary) {
			sE2E(3 * f0 + 2, -1);
			sE2E(3 * f3 + 0, -1);
		}
		else {
			const int e1p = E2E[dedge_prev_3(e1)],
				e1n = E2E[dedge_next_3(e1)];
			sE2E(3 * f0 + 2, 3 * f1 + 0);
			sE2E(3 * f1 + 1, e1n);
			sE2E(3 * f1 + 2, 3 * f2 + 0);
			sE2E(3 * f2 + 1, e1p);
			sE2E(3 * f2 + 2, 3 * f3 + 0);
		}
#undef sE2E

		/* Update V2E */
		V2E[v0] = 3 * f0 + 2;
		V2E[vn] = 3 * f0 + 0;
		V2E[v1] = 3 * f3 + 1;
		V2E[v0p] = 3 * f0 + 1;
		if (!is_boundary)
			V2E[v1p] = 3 * f1 + 2;

		auto schedule = [&](int f) {
			for (int i = 0; i<3; ++i) {
				double length = (V.col(F(i, f)) - V.col(F((i + 1) % 3, f))).squaredNorm();
				if (length > maxLength)
					queue.push(Edge(length, f * 3 + i));
			}
		};

		schedule(f0);
		if (!is_boundary) {
			schedule(f2);
			schedule(f1);
		};
		schedule(f3);
	}
	F.conservativeResize(F.rows(), nF);
	V.conservativeResize(V.rows(), nV);
	V2E.conservativeResize(nV);
	boundary.conservativeResize(nV);
	nonmanifold.conservativeResize(nV);
	E2E.conservativeResize(nF * 3);


}

#include "parametrizer.hpp"
#include "disajoint-tree.hpp"
#include "field-math.hpp"

void subdivide_diff(MatrixXi &F,
	MatrixXd &V, 
	MatrixXd &N,
	MatrixXd &Q,
	MatrixXd &O,
	VectorXi &V2E, VectorXi &E2E,
	VectorXi &boundary, VectorXi &nonmanifold, 
	std::vector<Vector2i>& edge_diff,
	std::vector<DEdge>& edge_values,
	std::vector<Vector3i>& face_edgeOrients,
	std::vector<Vector3i>& face_edgeIds,
	std::map<int, int>& singularities
	) {
	struct EdgeInfo
	{
		int e2eid;
		DEdge edge;
	};

	std::vector<EdgeInfo> queue;
	int qid = 0;
	bool debug = false;
	auto sanity = [&](int c) {
		printf("check sanity %d\n", c);
		for (int i = 0; i < face_edgeIds.size(); ++i) {
			for (int j = 0; j < 3; ++j) {
				if (edge_values[face_edgeIds[i][j]] != DEdge(F(j, i), F((j + 1) % 3, i))) {
					printf("edge id wrong!\n");
                    exit(0);
				}
			}
		}
		for (int i = 0; i < face_edgeOrients.size(); ++i) {
			if (singularities.count(i))
				continue;
			for (int j = 0; j < 3; ++j) {
				int od0, od1;
				int orient0 = face_edgeOrients[i][j];
				int orient1 = face_edgeOrients[i][(j + 1) % 3];
				int v0 = F(j, i), v1 = F((j + 1) % 3, i), v2 = F((j + 2) % 3, i);
				if (v0 < v1) {
					auto value1 = compat_orientation_extrinsic_index_4(Q.col(v0), N.col(v0), Q.col(F(0, i)), N.col(F(0, i)));
					od0 = (value1.second - value1.first + 4) % 4;
				}
				else {
					auto value1 = compat_orientation_extrinsic_index_4(Q.col(v1), N.col(v1), Q.col(F(0, i)), N.col(F(0, i)));
					od0 = (value1.second - value1.first + 6) % 4;
				}
				if (v1 < v2) {
					auto value1 = compat_orientation_extrinsic_index_4(Q.col(v1), N.col(v1), Q.col(F(0, i)), N.col(F(0, i)));
					od1 = (value1.second - value1.first + 4) % 4;
				}
				else {
					auto value1 = compat_orientation_extrinsic_index_4(Q.col(v2), N.col(v2), Q.col(F(0, i)), N.col(F(0, i)));
					od1 = (value1.second - value1.first + 6) % 4;
				}

				if ((orient0 - od0 + 4) % 4 != (orient1 - od1 + 4) % 4) {
					printf("%d %d %d %d\n", orient0, od0, orient1, od1);
					printf("orient wrong!\n");
                    exit(0);
                }
			}
		}
		for (int i = 0; i < face_edgeIds.size(); ++i) {
			Vector2i total_diff(0, 0);
			for (int j = 0; j < 3; ++j) {
				total_diff += rshift90(edge_diff[face_edgeIds[i][j]], face_edgeOrients[i][j]);
				if ((c || singularities.count(i)) &&
					(abs(edge_diff[face_edgeIds[i][j]][0]) > 1 || abs(edge_diff[face_edgeIds[i][j]][1]) > 1))
				{
					printf("long edge %d!\n", singularities.count(i));
					printf("%d %d\n", i, j);
                    exit(0);
				}
			}
			if (total_diff != Vector2i(0, 0)) {
				printf("non zero diff\n");
                exit(0);
			}
		}
		printf("finish check!\n");
	};

	auto AnalyzeOrients = [&](int f0, int orients[3]) {
		int vo = F(0, f0);
		for (int j = 0; j < 3; ++j) {
			int od0, od1;
			int v0 = F(j, f0), v1 = F((j + 1) % 3, f0);
			if (v0 < v1) {
				auto value1 = compat_orientation_extrinsic_index_4(Q.col(v0), N.col(v0), Q.col(vo), N.col(vo));
				od0 = (value1.second - value1.first + 4) % 4;
			}
			else {
				auto value1 = compat_orientation_extrinsic_index_4(Q.col(v1), N.col(v1), Q.col(vo), N.col(vo));
				od0 = (value1.second - value1.first + 6) % 4;
			}
			orients[j] = od0;
		}
	};

	sanity(0);
	for (int i = 0; i < E2E.size(); ++i) {
		if (singularities.count(i / 3) || singularities.count(E2E[i] / 3)) {
			continue;
		}
		int v0 = F(i % 3, i / 3), v1 = F((i + 1) % 3, i / 3);
		if (nonmanifold[v0] || nonmanifold[v1])
			continue;
		Vector2i diff = edge_diff[face_edgeIds[i/3][i%3]];
		if (abs(diff[0]) > 1 || abs(diff[1]) > 1) {
			int other = E2E[i];
			if (other == -1 || other > i) {
				EdgeInfo info;
				info.e2eid = i;
				info.edge = DEdge(v0, v1);
				queue.push_back(info);
			}
		}
	}

	int nV = V.cols(), nF = F.cols(), nSplit = 0;
	int counter = 0;

	while (qid != queue.size()) {
		counter += 1;
		EdgeInfo edge = queue[qid];
		qid++;
		int e0 = edge.e2eid, e1 = E2E[e0];
		if (edge.edge != edge_values[face_edgeIds[e0 / 3][e0 % 3]]) {
			continue;
		}
		bool is_boundary = e1 == -1;
		int f0 = e0 / 3, f1 = is_boundary ? -1 : (e1 / 3);
		int f2 = is_boundary ? -1 : (nF++);
		if (f2 != -1) {
			face_edgeOrients.push_back(Vector3i());
			face_edgeIds.push_back(Vector3i());
		}
		int f3 = nF++;
		face_edgeOrients.push_back(Vector3i());
		face_edgeIds.push_back(Vector3i());
		int v0 = F(e0 % 3, f0), v0p = F((e0 + 2) % 3, f0), v1 = F((e0 + 1) % 3, f0);
		int v1p = is_boundary ? -1 : F((e1 + 2) % 3, f1);

		int eid = face_edgeIds[f0][e0 % 3];
		int eid01 = face_edgeIds[f0][(e0 + 1) % 3];
		int eid02 = face_edgeIds[f0][(e0 + 2) % 3];
		int eid11 = face_edgeIds[f1][(e1 + 1) % 3];
		int eid12 = face_edgeIds[f1][(e1 + 2) % 3];

		int eid0, eid1, eid0p, eid1p;

		int vn = nV++;
		nSplit++;
		/* Update V */
		if (nV > V.cols()) {
			V.conservativeResize(V.rows(), V.cols() * 2);
			N.conservativeResize(N.rows(), N.cols() * 2);
			Q.conservativeResize(Q.rows(), Q.cols() * 2);
			O.conservativeResize(O.rows(), O.cols() * 2);
			V2E.conservativeResize(V.cols());
			boundary.conservativeResize(V.cols());
			nonmanifold.conservativeResize(V.cols());
		}

		/* Update V */
		V.col(vn) = (V.col(v0) + V.col(v1)) * 0.5;
		N.col(vn) = N.col(v0);
		auto value = compat_orientation_extrinsic_4(Q.col(v0), N.col(v0), Q.col(v1), N.col(v1));
		Q.col(vn) = Q.col(v0);
		O.col(vn) = (O.col(v0) + O.col(v1)) * 0.5;

		auto id1 = compat_orientation_extrinsic_index_4(Q.col(v0), N.col(v0), Q.col(vn), N.col(vn));
		auto id2 = compat_orientation_extrinsic_index_4(Q.col(vn), N.col(vn), Q.col(v1), N.col(v1));
		auto id3 = compat_orientation_extrinsic_index_4(Q.col(v0), N.col(v0), Q.col(v1), N.col(v1));
		int t1 = id1.second - id1.first;
		int t2 = id2.second - id2.first;
		if ((t1 + t2 + 8) % 4 != (id3.second - id3.first + 4) % 4) {
			printf("inconsistent orients!\n");
		}
		nonmanifold[vn] = false;
		boundary[vn] = is_boundary;
		int orient0, orient1, orient0p, orient1p;
		eid0 = eid;
		edge_values[eid0] = DEdge(v0, vn);
		
		eid1 = edge_values.size();
		edge_values.push_back(DEdge(vn, v1));
		auto index = compat_orientation_extrinsic_index_4(Q.col(v0), N.col(v0), Q.col(v1), N.col(v1));
		if (v0 < v1) {
			auto diff = edge_diff[eid];
			edge_diff[eid0] = diff / 2;
			edge_diff.push_back(rshift90(diff / 2 - diff, (index.second - index.first + 4) % 4));
		}
		else {
			auto diff = edge_diff[eid];
			edge_diff[eid0] = rshift90(-diff / 2, (index.first - index.second + 4) % 4);
			edge_diff.push_back(diff - diff / 2);
		}

		eid0p = edge_values.size();
		edge_values.push_back(DEdge(vn, v0p));

		/* Update F and E2E */
		int f20, f30;
		if (nF > F.cols()) {
			F.conservativeResize(F.rows(), std::max(nF, (int)F.cols() * 2));
			E2E.conservativeResize(F.cols() * 3);
		}

		/* Update F */
		F.col(f0) << vn, v0p, v0;
		int vo = F(0, f0);
		face_edgeIds[f0] = Vector3i(eid0p, eid02, eid0);

		int orients[3];
		AnalyzeOrients(f0, orients);
		face_edgeOrients[f0] = Vector3i(orients[0], orients[1], orients[2]);
		edge_diff.push_back(rshift90(-rshift90(edge_diff[eid02], orients[1]) - rshift90(edge_diff[eid0], orients[2]), (4 - orients[0]) % 4));

		if (!is_boundary) {
			F.col(f1) << vn, v0, v1p;
			AnalyzeOrients(f1, orients);
			eid1p = edge_values.size();
			edge_values.push_back(DEdge(vn, v1p));

			face_edgeOrients[f1] = Vector3i(orients[0], orients[1], orients[2]);
			face_edgeIds[f1] = (Vector3i(eid0, eid11, eid1p));
			edge_diff.push_back(rshift90(-rshift90(edge_diff[eid0], orients[0]) - rshift90(edge_diff[eid11], orients[1]), (4 - orients[2]) % 4));

			F.col(f2) << vn, v1p, v1;
			AnalyzeOrients(f2, orients);
			face_edgeIds[f2] = (Vector3i(eid1p, eid12, eid1));
			face_edgeOrients[f2] = Vector3i(orients[0], orients[1], orients[2]);
		}
		F.col(f3) << vn, v1, v0p;
		AnalyzeOrients(f3, orients);
		face_edgeIds[f3] = (Vector3i(eid1, eid01, eid0p));
		face_edgeOrients[f3] = Vector3i(orients[0], orients[1], orients[2]);

		/* Update E2E */
		const int e0p = E2E[dedge_prev_3(e0)],
			e0n = E2E[dedge_next_3(e0)];
		const int e1p = E2E[dedge_prev_3(e1)],
			e1n = E2E[dedge_next_3(e1)];

#define sE2E(a, b) E2E[a] = b; if (b != -1) E2E[b] = a;
		sE2E(3 * f0 + 0, 3 * f3 + 2);
		sE2E(3 * f0 + 1, e0p);
		sE2E(3 * f3 + 1, e0n);
		if (is_boundary) {
			sE2E(3 * f0 + 2, -1);
			sE2E(3 * f3 + 0, -1);
		}
		else {
			sE2E(3 * f0 + 2, 3 * f1 + 0);
			sE2E(3 * f1 + 1, e1n);
			sE2E(3 * f1 + 2, 3 * f2 + 0);
			sE2E(3 * f2 + 1, e1p);
			sE2E(3 * f2 + 2, 3 * f3 + 0);
		}
#undef sE2E

		/* Update V2E */
		V2E[v0] = 3 * f0 + 2;
		V2E[vn] = 3 * f0 + 0;
		V2E[v1] = 3 * f3 + 1;
		V2E[v0p] = 3 * f0 + 1;
		if (!is_boundary)
			V2E[v1p] = 3 * f1 + 2;

		auto schedule = [&](int f) {
			for (int i = 0; i<3; ++i) {
				Vector2i diff = edge_diff[face_edgeIds[f][i]];
				if (abs(diff[0]) > 1 || abs(diff[1]) > 1) {
					EdgeInfo info;
					info.e2eid = f * 3 + i;
					info.edge = edge_values[face_edgeIds[f][i]];
					queue.push_back(info);

				}
			}
		};

		schedule(f0);
		if (!is_boundary) {
			schedule(f2);
			schedule(f1);
		};
		schedule(f3);
	}
	sanity(1);

	F.conservativeResize(F.rows(), nF);
	V.conservativeResize(V.rows(), nV);
	N.conservativeResize(N.rows(), nV);
	Q.conservativeResize(Q.rows(), nV);
	O.conservativeResize(O.rows(), nV);

	V2E.conservativeResize(nV);
	boundary.conservativeResize(nV);
	nonmanifold.conservativeResize(nV);
	E2E.conservativeResize(nF * 3);


}
