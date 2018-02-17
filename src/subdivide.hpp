#include <Eigen/Core>
#include <Eigen/Dense>
#include "parametrizer.hpp"
using namespace Eigen;

void subdivide(MatrixXi &F, MatrixXd &V, VectorXi &V2E, VectorXi &E2E, VectorXi &boundary,
               VectorXi &nonmanifold, double maxLength);

void subdivide_diff(MatrixXi &F, MatrixXd &V, MatrixXd &N, MatrixXd &Q, MatrixXd &O, VectorXi &V2E,
                    VectorXi &E2E, VectorXi &boundary, VectorXi &nonmanifold,
                    std::vector<Vector2i> &edge_diff, std::vector<DEdge> &edge_values,
                    std::vector<Vector3i> &face_edgeOrients, std::vector<Vector3i> &face_edgeIds,
                    std::map<int, int> &singularities, int max_len);
