#include <Eigen/Core>
#include <Eigen/Dense>

using namespace Eigen;

void subdivide(MatrixXi &F, MatrixXf &V, VectorXi &V2E, VectorXi &E2E,
	VectorXi &boundary, VectorXi &nonmanifold, float maxLength);
