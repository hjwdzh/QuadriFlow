#ifndef CURVATURE_H_
#define CURVATURE_H_

#include <Eigen/Core>
#include <Eigen/Dense>

void ComputeCurvature(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::VectorXd& K);

#endif
