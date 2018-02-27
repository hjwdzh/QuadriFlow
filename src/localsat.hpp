#ifndef __LOCAL_SAT_H
#define __LOCAL_SAT_H

#include <Eigen/Core>
#include <vector>

using namespace Eigen;
void ExportLocalSat(std::vector<Vector2i> &edge_diff, const std::vector<Vector3i> &face_edgeIds,
                    const std::vector<Vector3i> &face_edgeOrients, const MatrixXi &F,
                    const VectorXi &V2E, const VectorXi &E2E);

#endif
