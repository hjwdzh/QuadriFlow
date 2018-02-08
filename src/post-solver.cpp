//
//  post-solver.cpp
//  parametrize
//
//  Created by Jingwei on 2/5/18.
//

#include "post-solver.hpp"
#include "ceres/ceres.h"
#include "ceres/rotation.h"

template <typename T, typename T2>
T Dot(const T a[3], const T2 b[3]) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

template <typename T>
T Length2(const T a[3]) {
    return Dot(a, a);
}

struct FaceConstraint {
    FaceConstraint(double coeff_length, double coeff_bias, Vector3d normal[4], double bias[4],
                   double length)
        : coeff_length(coeff_length),
          coeff_bias(coeff_bias),
          length0(length),
          bias0{
              bias[0],
              bias[1],
              bias[2],
              bias[3],
          },
          normal0{
              normal[0],
              normal[1],
              normal[2],
              normal[3],
          } {}

    template <typename T>
    bool operator()(const T* p0, T* r) const {
        const T* p[] = {p0, p0 + 3, p0 + 6, p0 + 9};
        for (int k = 0; k < 4; ++k) {
            auto pc = p[k];
            auto pa = p[(k + 1) % 4];
            auto pb = p[(k + 3) % 4];

            T a[3]{pa[0] - pc[0], pa[1] - pc[1], pa[2] - pc[2]};
            T b[3]{pb[0] - pc[0], pb[1] - pc[1], pb[2] - pc[2]};
            r[3 * k + 0] = coeff_length * (ceres::sqrt(Length2(a)) - length0);

            T normal[3];
            ceres::CrossProduct(a, b, normal);
            T l2normal = ceres::sqrt(Length2(normal));
            if (l2normal == T()) continue;
            for (int i = 0; i < 3; ++i) normal[i] /= l2normal;
            r[3 * k + 1] = ceres::acos(Dot(normal, &normal0[k][0]));
            r[3 * k + 2] = coeff_bias * (Dot(pc, normal) - bias0[k]);
        }
        return true;
    }

    static ceres::CostFunction* create(double alpha, double beta, Vector3d normal[4],
                                       double bias[4], double length) {
        return new ceres::AutoDiffCostFunction<FaceConstraint, 3 * 4, 4 * 3>(
            new FaceConstraint(alpha, beta, normal, bias, length));
    }

    double coeff_length;
    double coeff_bias;
    double length0;

    double bias0[4];
    Vector3d normal0[4];
};

void optimize_quad_positions(std::vector<Vector3d>& O_quad, std::vector<Vector3d>& N_quad,
                             std::vector<Vector3d>& Q_quad, std::vector<Vector4i>& F_quad,
                             VectorXi& V2E_quad, VectorXi& E2E_quad, MatrixXd& V, MatrixXd& N,
                             MatrixXd& Q, MatrixXd& O, MatrixXi& F, VectorXi& V2E, VectorXi& E2E,
                             DisajointTree& disajoint_tree, double reference_length) {
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
    printf("Reference length: %.2f\n", reference_length);

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

    int n_quad = Q_quad.size();
    int n_trig = Q.size();
    std::vector<double> B_quad(n_quad);
    std::vector<int> B_weight(n_quad);
    for (int vtrig = 0; vtrig < n_trig; ++vtrig) {
        int vquad = disajoint_tree.Index(vtrig);
        double b = N.col(vtrig).dot(O.col(vtrig));
        B_quad[vquad] += b;
        B_weight[vquad] += 1;
    }
    for (int vquad = 0; vquad < n_quad; ++vquad) B_quad[vquad] /= B_weight[vquad];

    ceres::Problem problem;
    std::vector<double> solution(n_quad * 3);
    for (int vquad = 0; vquad < n_quad; ++vquad) {
        solution[3 * vquad + 0] = O_quad[vquad][0];
        solution[3 * vquad + 1] = O_quad[vquad][1];
        solution[3 * vquad + 2] = O_quad[vquad][2];
    }

    for (int fquad = 0; fquad < F_quad.size(); ++fquad) {
        auto v = F_quad[fquad];

        double bias[4];
        Vector3d normal[4];
        std::vector<double*> var(12);
        for (int k = 0; k < 4; ++k) {
            var[3 * k + 0] = &solution[3 * v[k] + 0];
            var[3 * k + 1] = &solution[3 * v[k] + 1];
            var[3 * k + 2] = &solution[3 * v[k] + 2];
            bias[k] = B_quad[k];
            normal[k] = N_quad[k];
        }
        ceres::CostFunction* cost_function =
            FaceConstraint::create(0.01, 0.5, normal, bias, reference_length);
        problem.AddResidualBlock(cost_function, nullptr, var);
    }

    ceres::Solver::Options options;
    options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
    options.minimizer_progress_to_stdout = true;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);

    std::cout << summary.BriefReport() << std::endl;

    return;

    for (int vquad = 0; vquad < n_quad; ++vquad) {
        O_quad[vquad][0] = solution[3 * vquad + 0];
        O_quad[vquad][1] = solution[3 * vquad + 1];
        O_quad[vquad][2] = solution[3 * vquad + 2];
    }
}
