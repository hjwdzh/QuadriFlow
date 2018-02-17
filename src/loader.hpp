#ifndef __LOADER_H
#define __LOADER_H

#include <Eigen/Core>
#include <vector>

using namespace Eigen;

void load(const char* filename, MatrixXd& V, MatrixXi& F);

#endif
