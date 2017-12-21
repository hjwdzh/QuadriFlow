#include <Eigen/Core>
#include <Eigen/Dense>

using namespace Eigen;

void subdivide(MatrixXi &F, MatrixXd &V, VectorXi &V2E, VectorXi &E2E,
	VectorXi &boundary, VectorXi &nonmanifold, double maxLength);
