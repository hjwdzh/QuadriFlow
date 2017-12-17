#include "MergeVertex.h"

#include "CompareKey.h"

#include <map>
#include <vector>

void merge_close(MatrixXf& V, MatrixXi& F, float threshold)
{
	std::map<Key3f, int> vid_maps;
	std::vector<int> vid_compress(V.cols());
	for (int i = 0; i < V.cols(); ++i) {
		Key3f key(V(0, i), V(1, i), V(2, i), threshold);
		if (vid_maps.count(key)) {
			vid_compress[i] = vid_maps[key];
		}
		else {
			V.col(vid_maps.size()) = V.col(i);
			vid_compress[i] = vid_maps.size();
			vid_maps[key] = vid_compress[i];
		}
	}
	printf("Compress from %d to %d...\n", V.cols(), vid_maps.size());
	MatrixXf newV(3, vid_maps.size());
	memcpy(newV.data(), V.data(), sizeof(float) * 3 * vid_maps.size());
	V = std::move(newV);
	for (int j = 0; j < 3; ++j) {
		for (int i = 0; i < F.cols(); ++i) {
			F(j, i) = vid_compress[F(j, i)];
		}
	}
}