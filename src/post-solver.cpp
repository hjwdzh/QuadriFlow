//
//  post-solver.cpp
//  parametrize
//
//  Created by Jingwei on 2/5/18.
//

#include "post-solver.hpp"
#include "ceres/ceres.h"
#include "ceres/rotation.h"

struct FaceConstraint {
    FaceConstraint(double alpha, double beta, Vector3d normal[4], double bias, double length)
        : alpha(alpha),
          beta(beta),
          normal0{normal[0], normal[1], normal[2], normal[3]},
          bias0(bias),
          length0(length) {}

    template <typename T>
    T dot(const T a[3], const T b[3]) {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    }

    template <typename T>
    T length2(const T a[3]) {
        return dot(a, a);
    }

    template <typename T>
    bool operator()(const T* p0, T* residuals) const {
        const T* p[] = {p0, p0 + 3, p0 + 6, p0 + 9};
        for (int k = 0; k < 4; ++k) {
            auto pc = p[k];
            auto pa = p[(k + 1) % 4];
            auto pb = p[(k + 2) % 4];

            T a[3]{pa[0] - pc[0], pa[1] - pc[1], pa[2] - pc[2]};
            T b[3]{pb[0] - pc[0], pb[1] - pc[1], pb[2] - pc[2]};
            T normal[3];
            ceres::CrossProduct(a, b, normal);
            // length2(normal);
            // T l2normal = ceres::sqrt();
        }
        return true;
    }

    static ceres::CostFunction* create(double alpha, double beta, Vector3d normal[4], double bias,
                                       double length) {
        return new ceres::AutoDiffCostFunction<FaceConstraint, 1, 4 * 3>(
            new FaceConstraint(alpha, beta, normal, bias, length));
    }

    double alpha;
    double beta;
    Vector3d normal0[4];
    double bias0;
    double length0;
};

void optimize_quad_positions(std::vector<Vector3d>& O_quad, std::vector<Vector3d>& N_quad,
                             std::vector<Vector3d>& Q_quad, std::vector<Vector4i>& F_quad,
                             VectorXi& V2E_quad, VectorXi& E2E_quad, MatrixXd& V, MatrixXd& N,
                             MatrixXd& Q, MatrixXd& O, MatrixXi& F, VectorXi& V2E, VectorXi& E2E,
                             DisajointTree& disajoint_tree) {
    // Information for the quad mesh
    printf("Quad mesh info:\n");
    printf("Number of vertices with normals and orientations: %d = %d = %d\n", (int)O_quad.size(),
           (int)N_quad.size(), (int)Q_quad.size());
    printf("Number of faces: %d\n", (int)F_quad.size());
    printf("Number of directed edges: %d\n", (int)E2E_quad.size());
    // Information for the original mesh
    printf("Triangle mesh info:\n");
    printf(
        "Number of vertices with normals,"
        "orientations and associated quad positions: "
        "%d = %d = %d = %d\n",
        (int)V.cols(), (int)N.cols(), (int)Q.cols(), (int)O.cols());
    printf("Number of faces: %d\n", (int)F.cols());
    printf("Number of directed edges: %d\n", (int)E2E.size());

    /* initial quad flips
     currently there are many flips, (82 flips in hand.obj)
     By uncommenting the ComputePosition() function call in Parametrizer.cpp, (current post linear
     solver) the flip number will be reduced to 14.

     As we discussed, we hope to implement this function to replace the previous solver
     */
    int flip_count = 0;
    for (int i = 0; i < F_quad.size(); ++i) {
        bool flipped = false;
        for (int j = 0; j < 4; ++j) {
            int v1 = F_quad[i][j];
            int v2 = F_quad[i][(j + 1) % 4];
            int v3 = F_quad[i][(j + 3) % 4];

            Vector3d face_norm = (O_quad[v2] - O_quad[v1]).cross(O_quad[v3] - O_quad[v1]);
            Vector3d vertex_norm = N_quad[v1];
            if (face_norm.dot(vertex_norm) < 0) {
                flipped = true;
            }
        }
        if (flipped) {
            flip_count++;
        }
    }
    printf("Flipped Quads: %d\n", flip_count);
}
