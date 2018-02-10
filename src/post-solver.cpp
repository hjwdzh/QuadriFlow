//
//  post-solver.cpp
//  parametrize
//
//  Created by Jingwei on 2/5/18.
//
#include <cmath>
#include <cstdio>
#include <string>
#include "ceres/ceres.h"
#include "ceres/rotation.h"

#include "post-solver.hpp"
#include "serialize.hpp"

const double COEFF_AREA = 1;
const double COEFF_BIAS = 1;
const double COEFF_NORMAL = 1;
const int N_ITER = 100;

template <typename T, typename T2>
T DotProduct(const T a[3], const T2 b[3]) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

template <typename T>
T Length2(const T a[3]) {
    return DotProduct(a, a);
}

bool DEBUG = 0;
struct FaceConstraint {
    FaceConstraint(double coeff_area, double coeff_bias, double coeff_normal, Vector3d normal[4],
                   double bias[4], double length)
        : coeff_area(coeff_area),
          coeff_bias(coeff_bias),
          coeff_normal(coeff_normal),
          area0(length * length),
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
    bool operator()(const T* p0, const T* p1, const T* p2, const T* p3, T* r) const {
        const T* p[] = {p0, p1, p2, p3};
        r[8] = T();
        for (int k = 0; k < 4; ++k) {
            auto pc = p[k];
            auto pa = p[(k + 1) % 4];
            auto pb = p[(k + 3) % 4];

            T a[3]{pa[0] - pc[0], pa[1] - pc[1], pa[2] - pc[2]};
            T b[3]{pb[0] - pc[0], pb[1] - pc[1], pb[2] - pc[2]};

            T normal[3];
            ceres::CrossProduct(a, b, normal);
            T area = ceres::sqrt(Length2(normal));
            r[8] += area;

            if (area == T()) continue;
            for (int i = 0; i < 3; ++i) normal[i] /= area;
            T degree = ceres::acos(DotProduct(normal, &normal0[k][0]));
            r[2 * k + 0] = coeff_normal * degree * degree;
            r[2 * k + 1] = coeff_bias * (DotProduct(pc, normal) - bias0[k]);
        }
        r[8] = coeff_area * (r[8] / (4.0 * area0) - 1.0);
        return true;
    }

    static ceres::CostFunction* create(double coeff_area, double coeff_bias, double coeff_normal,
                                       Vector3d normal[4], double bias[4], double length) {
        return new ceres::AutoDiffCostFunction<FaceConstraint, 9, 3, 3, 3, 3>(
            new FaceConstraint(coeff_area, coeff_bias, coeff_normal, normal, bias, length));
    }

    double coeff_area;
    double coeff_bias;
    double coeff_normal;
    double area0;

    double bias0[4];
    Vector3d normal0[4];
};

void solve(std::vector<Vector3d>& O_quad, std::vector<Vector3d>& N_quad,
           std::vector<Vector3d>& Q_quad, std::vector<Vector4i>& F_quad,
           std::vector<double>& B_quad, MatrixXd& V, MatrixXd& N, MatrixXd& Q, MatrixXd& O,
           MatrixXi& F, double reference_length) {
    int n_quad = Q_quad.size();

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
        for (int k = 0; k < 4; ++k) {
            bias[k] = B_quad[v[k]];
            normal[k] = N_quad[v[k]];
        }
        ceres::CostFunction* cost_function = FaceConstraint::create(
            COEFF_AREA, COEFF_BIAS, COEFF_NORMAL, normal, bias, reference_length);
        /*
        double r[9];
        DEBUG = fquad < 2;
        FaceConstraint(1, 1, normal, bias, reference_length)(
            &solution[3 * v[0]], &solution[3 * v[1]], &solution[3 * v[2]], &solution[3 * v[3]], r);
        for (int i = 0; i < 9; ++i) init_error += r[i] * r[i];
        */
        problem.AddResidualBlock(cost_function, nullptr, &solution[3 * v[0]], &solution[3 * v[1]],
                                 &solution[3 * v[2]], &solution[3 * v[3]]);
    }

    ceres::Solver::Options options;
    options.num_threads = 1;
    options.max_num_iterations = N_ITER;
    options.initial_trust_region_radius = 1;
    options.linear_solver_type = ceres::CGNR;
    options.minimizer_progress_to_stdout = true;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);

    std::cout << summary.BriefReport() << std::endl;

    for (int vquad = 0; vquad < n_quad; ++vquad) {
        O_quad[vquad][0] = solution[3 * vquad + 0];
        O_quad[vquad][1] = solution[3 * vquad + 1];
        O_quad[vquad][2] = solution[3 * vquad + 2];
    }

    return;
}

void optimize_quad_positions(std::vector<Vector3d>& O_quad, std::vector<Vector3d>& N_quad,
                             std::vector<Vector3d>& Q_quad, std::vector<Vector4i>& F_quad,
                             VectorXi& V2E_quad, std::vector<int>& E2E_quad, MatrixXd& V, MatrixXd& N,
                             MatrixXd& Q, MatrixXd& O, MatrixXi& F, VectorXi& V2E, VectorXi& E2E,
                             DisajointTree& disajoint_tree, double reference_length,
                             bool just_serialize) {
    printf("Quad mesh info:\n");
    printf("Number of vertices with normals and orientations: %d = %d = %d\n", (int)O_quad.size(),
           (int)N_quad.size(), (int)Q_quad.size());
    printf("Number of faces: %d\n", (int)F_quad.size());
    printf("Number of directed edges: %d\n", (int)E2E_quad.size());
    // Information for the original mesh
    printf("Triangle mesh info:\n");
    printf(
        "Number of vertices with normals, "
        "orientations and associated quad positions: "
        "%d = %d = %d = %d\n",
        (int)V.cols(), (int)N.cols(), (int)Q.cols(), (int)O.cols());
    printf("Number of faces: %d\n", (int)F.cols());
    printf("Number of directed edges: %d\n", (int)E2E.size());
    printf("Reference length: %.2f\n", reference_length);

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

    int n_quad = O_quad.size();
    int n_trig = O.cols();
    std::vector<double> B_quad(n_quad);  // Average bias for quad vertex
    std::vector<int> B_weight(n_quad);

    printf("ntrig: %d, disjoint_tree.size: %d\n", n_trig, (int)disajoint_tree.indices.size());
    for (int vtrig = 0; vtrig < n_trig; ++vtrig) {
        int vquad = disajoint_tree.Index(vtrig);
        double b = N_quad[vquad].dot(O.col(vtrig));
        B_quad[vquad] += b;
        B_weight[vquad] += 1;
    }
    for (int vquad = 0; vquad < n_quad; ++vquad) {
        assert(B_weight[vquad]);
        B_quad[vquad] /= B_weight[vquad];
    }

    printf("just_serialize: %d\n", just_serialize);
    if (just_serialize) {
        puts("Save parameters to post.bin for optimization");
        FILE* out = fopen("post.bin", "wb");
        assert(out);
        Save(out, O_quad);
        Save(out, N_quad);
        Save(out, Q_quad);
        Save(out, F_quad);
        Save(out, B_quad);
        Save(out, V);
        Save(out, N);
        Save(out, Q);
        Save(out, O);
        Save(out, F);
        Save(out, reference_length);
        fclose(out);
    } else {
        puts("Start post optimization");
        solve(O_quad, N_quad, Q_quad, F_quad, B_quad, V, N, Q, O, F, reference_length);
    }
}

#ifdef POST_SOLVER

void SaveObj(const std::string& fname, std::vector<Vector3d> O_quad,
             std::vector<Vector4i> F_quad) {
    std::ofstream os(fname);
    for (int i = 0; i < (int)O_quad.size(); ++i) {
        os << "v " << O_quad[i][0] << " " << O_quad[i][1] << " " << O_quad[i][2] << "\n";
    }
    for (int i = 0; i < (int)F_quad.size(); ++i) {
        os << "f " << F_quad[i][0] + 1 << " " << F_quad[i][1] + 1 << " " << F_quad[i][2] + 1 << " "
           << F_quad[i][3] + 1 << "\n";
    }
    os.close();
}

int main() {
    std::vector<Vector3d> O_quad;
    std::vector<Vector3d> N_quad;
    std::vector<Vector3d> Q_quad;
    std::vector<Vector4i> F_quad;
    std::vector<double> B_quad;
    MatrixXd V;
    MatrixXd N;
    MatrixXd Q;
    MatrixXd O;
    MatrixXi F;
    double reference_length;

    puts("Read parameters from post.bin");
    FILE* in = fopen("post.bin", "rb");
    assert(in);
    Read(in, O_quad);
    Read(in, N_quad);
    Read(in, Q_quad);
    Read(in, F_quad);
    Read(in, B_quad);
    Read(in, V);
    Read(in, N);
    Read(in, Q);
    Read(in, O);
    Read(in, F);
    Read(in, reference_length);
    fclose(in);
    printf("reference_length: %.2f\n", reference_length);

    int n_flip = 0;
    double sum_degree = 0;
    for (int i = 0; i < F_quad.size(); ++i) {
        bool flipped = false;
        for (int j = 0; j < 4; ++j) {
            int v1 = F_quad[i][j];
            int v2 = F_quad[i][(j + 1) % 4];
            int v3 = F_quad[i][(j + 3) % 4];

            Vector3d face_norm =
                (O_quad[v2] - O_quad[v1]).cross(O_quad[v3] - O_quad[v1]).normalized();
            Vector3d vertex_norm = N_quad[v1];
            if (face_norm.dot(vertex_norm) < 0) {
                flipped = true;
            }
            double degree = std::acos(face_norm.dot(vertex_norm));
            assert(degree >= 0);
            // printf("cos theta = %.2f\n", degree);
            sum_degree += degree * degree;
        }
        n_flip += flipped;
    }
    printf("n_flip: %d\nsum_degree: %.3f\n", n_flip, sum_degree);

    puts("Start post optimization");
    solve(O_quad, N_quad, Q_quad, F_quad, B_quad, V, N, Q, O, F, reference_length);
    SaveObj("postsolver.obj", O_quad, F_quad);

    n_flip = 0;
    sum_degree = 0;
    for (int i = 0; i < F_quad.size(); ++i) {
        bool flipped = false;
        for (int j = 0; j < 4; ++j) {
            int v1 = F_quad[i][j];
            int v2 = F_quad[i][(j + 1) % 4];
            int v3 = F_quad[i][(j + 3) % 4];

            Vector3d face_norm =
                (O_quad[v2] - O_quad[v1]).cross(O_quad[v3] - O_quad[v1]).normalized();
            Vector3d vertex_norm = N_quad[v1];
            if (face_norm.dot(vertex_norm) < 0) {
                flipped = true;
            }
            double degree = std::acos(face_norm.dot(vertex_norm));
            assert(degree >= 0);
            sum_degree += degree * degree;
        }
        n_flip += flipped;
    }
    printf("n_flip: %d\nsum_degree: %.3f\n", n_flip, sum_degree);
    return 0;
}

#endif
