//
//  post-solver.hpp
//  Parametrize
//
//  Created by Jingwei on 2/5/18.
//

#ifndef post_solver_h
#define post_solver_h

#include "disajoint-tree.hpp"
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
using namespace Eigen;


void optimize_quad_positions(std::vector<Vector3d>& O_quad,
                             std::vector<Vector3d>& N_quad,
                             std::vector<Vector3d>& Q_quad,
                             std::vector<Vector4i>& F_quad,
                             VectorXi& V2E_quad,
                             VectorXi& E2E_quad,
                             MatrixXd& V,
                             MatrixXd& N,
                             MatrixXd& Q,
                             MatrixXd& O,
                             MatrixXi& F,
                             VectorXi& V2E,
                             VectorXi& E2E,
                             DisajointTree& disajoint_tree);

#endif /* post_solver_h */
