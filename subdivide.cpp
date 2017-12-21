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
