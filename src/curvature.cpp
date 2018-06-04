#include "curvature.hpp"

/*
#include <igl/avg_edge_length.h>
#include <igl/cotmatrix.h>
#include <igl/invert_diag.h>
#include <igl/massmatrix.h>
#include <igl/parula.h>
#include <igl/per_corner_normals.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/principal_curvature.h>
using namespace Eigen;

void ComputeCurvature(const Eigen::MatrixXd& Vt, const Eigen::MatrixXi& Ft, Eigen::VectorXd& K) {
    // Compute integral of Gaussian curvature
    Eigen::MatrixXd V = Vt.transpose();
    Eigen::MatrixXi F = Ft.transpose();
    // Alternative discrete mean curvature
    MatrixXd HN;
    SparseMatrix<double> L, M, Minv;
    igl::cotmatrix(V, F, L);
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, M);
    igl::invert_diag(M, Minv);
    // Laplace-Beltrami of position
    HN = -Minv * (L * V);
    // Extract magnitude as mean curvature
    VectorXd H = HN.rowwise().norm();

    // Compute curvature directions via quadric fitting
    MatrixXd PD1, PD2;
    VectorXd PV1, PV2;
    igl::principal_curvature(V, F, PD1, PD2, PV1, PV2);
    K = Eigen::VectorXd(V.rows());
    // double min_rho = 1e30, max_rho = 0;
    for (int i = 0; i < PV1.size(); ++i) {
        double max_curv = std::max(std::abs(PV1[i]), std::abs(PV2[i]));
        double rho = 1.0 / max_curv;
        K[i] = rho;
    }
}

*/
