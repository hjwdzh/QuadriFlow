#include "field-math.hpp"

#include <Eigen/Core>
#include <vector>

using namespace Eigen;

void ExportLocalSat(std::vector<Vector2i> &edge_diff, std::vector<Vector3i> face_edgeIds,
                    std::vector<Vector3i> face_edgeOrients) {
    int flip_count = 0;
    for (int i = 0; i < face_edgeIds.size(); ++i) {
        Vector2i diff[3];
        for (int j = 0; j < 3; ++j) {
            diff[j] = rshift90(edge_diff[face_edgeIds[i][j]], face_edgeOrients[i][j]);
        }
        if (diff[0] + diff[1] + diff[2] != Vector2i::Zero()) {
            printf("Non zero!\n");
        }
        if (diff[0][0] * diff[1][1] - diff[0][1] * diff[1][0] < 0) {
            flip_count += 1;
        }
    }
}
