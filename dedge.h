#ifndef DEDGE_H_
#define DEDGE_H_

#include <Eigen/Core>
#include <Eigen/Dense>


using namespace Eigen;

inline int dedge_prev_3(int e) { return (e % 3 == 0) ? e + 2 : e - 1; }
inline int dedge_next_3(int e) { return (e % 3 == 2) ? e - 2 : e + 1; }

void compute_direct_graph(MatrixXf& V, MatrixXi& F, VectorXi& V2E,
	VectorXi& E2E, VectorXi& boundary, VectorXi& nonManifold);

#endif