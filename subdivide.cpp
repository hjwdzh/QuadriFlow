#include "subdivide.h"
#include "dedge.h"
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

#include "Parametrizer.h"
#include "DisajointTree.h"
#include "field_math.h"

void subdivide_diff(MatrixXi &F,
	MatrixXd &V, 
	MatrixXd &N,
	MatrixXd &Q,
	MatrixXd &O,
	VectorXi &V2E, VectorXi &E2E,
	VectorXi &boundary, VectorXi &nonmanifold, 
	std::vector<Vector2i>& edge_diff,
	std::vector<DEdge>& edge_values,
	std::map<DEdge,int>& edge_ids,
	std::vector<Vector4i>& edge_to_constraints,
	DisajointOrientTree& disajoint_orient_tree
	) {
	struct EdgeInfo
	{
		int eid;
		int e2eid;
		Vector2i diff;
	};

	std::queue<EdgeInfo> queue;

	for (int i = 0; i < E2E.size(); ++i) {
		int v0 = F(i % 3, i / 3), v1 = F((i + 1) % 3, i / 3);
		if (nonmanifold[v0] || nonmanifold[v1])
			continue;
		int eid = edge_ids[DEdge(v0, v1)];
		Vector2i diff = edge_diff[edge_ids[DEdge(v0, v1)]];
		if (abs(diff[0]) > 1 || abs(diff[1]) > 1) {
			int other = E2E[i];
			if (other == -1 || other > i) {
				EdgeInfo info;
				info.eid = eid;
				info.e2eid = i;
				info.diff = diff;
				queue.push(info);
			}
		}
	}

	int nV = V.cols(), nF = F.cols(), nSplit = 0;
	int counter = 0;
	while (!queue.empty()) {
		counter += 1;
		EdgeInfo edge = queue.front();
		queue.pop();
		int eid = edge.eid;
		int e0 = edge.e2eid, e1 = E2E[e0];
		bool is_boundary = e1 == -1;
		int f0 = e0 / 3, f1 = is_boundary ? -1 : (e1 / 3);
		int v0 = F(e0 % 3, f0), v0p = F((e0 + 2) % 3, f0), v1 = F((e0 + 1) % 3, f0);
		int v1p = is_boundary ? -1 : F((e1 + 2) % 3, f1);

		int eid01 = edge_ids[DEdge(v1, v0p)];
		int eid02 = edge_ids[DEdge(v0p, v0)];
		int eid11 = edge_ids[DEdge(v1, v1p)];
		int eid12 = edge_ids[DEdge(v1p, v0)];
		int orient0 = edge_to_constraints[eid][1];
		int orient1 = edge_to_constraints[eid][3];

		if (f0 == edge_to_constraints[eid][2])
			std::swap(orient0, orient1);

		DEdge eid0, eid1, eid0p, eid1p;
		int f_orient0 = disajoint_orient_tree.Orient(f0);
		int f_orient1 = disajoint_orient_tree.Orient(f1);

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
		N.col(vn) = (N.col(v0) + N.col(v1)).normalized();
		auto value = compat_orientation_extrinsic_4(Q.col(v0), N.col(v0), Q.col(v1), N.col(v1));
		Q.col(vn) = (value.first + value.second).normalized();
		O.col(vn) = (O.col(v0) + O.col(v1)) * 0.5;
		nonmanifold[vn] = false;
		boundary[vn] = is_boundary;
		
		eid0 = DEdge(vn, v0);
		edge_ids[eid0] = edge_values.size();
		edge_values.push_back(eid0);
		edge_diff.push_back(edge_diff[eid] / 2);

		eid1 = DEdge(vn, v1);
		edge_ids[eid1] = edge_values.size();
		edge_values.push_back(eid1);
		edge_diff.push_back(edge_diff[eid] / 2 - edge_diff[eid]);

		eid0p = DEdge(vn, v0p);
		edge_ids[eid0p] = edge_values.size();
		edge_values.push_back(eid0p);

		/* Update F and E2E */
		int f2 = is_boundary ? -1 : (nF++);
		int f3 = nF++;
		int f20, f30;
		if (nF > F.cols()) {
			F.conservativeResize(F.rows(), max(nF, (int)F.cols() * 2));
			E2E.conservativeResize(F.cols() * 3);
		}

		/* Update F */
		F.col(f0) << vn, v0p, v0;
		if (!is_boundary) {
			F.col(f1) << vn, v0, v1p;
			F.col(f2) << vn, v1p, v1;
			
			eid1p = DEdge(vn, v1p);
			edge_ids[eid1p] = edge_values.size();
			edge_values.push_back(eid1p);

			disajoint_orient_tree.parent.push_back(std::make_pair(disajoint_orient_tree.parent.size(), f_orient1));
		}
		F.col(f3) << vn, v1, v0p;
		disajoint_orient_tree.parent.push_back(std::make_pair(disajoint_orient_tree.parent.size(), f_orient0));

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
//				if (length > maxLength)
//					queue.push(Edge(length, f * 3 + i));
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
	N.conservativeResize(N.rows(), nV);
	Q.conservativeResize(Q.rows(), nV);
	O.conservativeResize(O.rows(), nV);

	V2E.conservativeResize(nV);
	boundary.conservativeResize(nV);
	nonmanifold.conservativeResize(nV);
	E2E.conservativeResize(nF * 3);


}
