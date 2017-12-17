#ifndef MERGE_VERTEX_H_
#define MERGE_VERTEX_H_

#include <Eigen/Core>
using namespace Eigen;

void merge_close(MatrixXf& V, MatrixXi& F, float threshold);

#endif