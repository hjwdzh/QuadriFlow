#include <glm/glm.hpp>
#include <cuda_runtime.h>
#include "AdjacentMatrix.h"

__device__ __host__ double cudaSignum(double value) {
	return std::copysign((double)1, value);
}

__device__ __host__ void
compat_orientation_extrinsic_4(const glm::dvec3 &q0, const glm::dvec3 &n0,
const glm::dvec3 &q1, const glm::dvec3 &n1, glm::dvec3& value1, glm::dvec3& value2) {
	const glm::dvec3 A[2] = { q0, glm::cross(n0, q0) };
	const glm::dvec3 B[2] = { q1, glm::cross(n1, q1) };

	double best_score = -1e10;
	int best_a = 0, best_b = 0;

	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < 2; ++j) {
			double score = std::abs(glm::dot(A[i], B[j]));
			if (score > best_score + 1e-6) {
				best_a = i;
				best_b = j;
				best_score = score;
			}
		}
	}
	const double dp = glm::dot(A[best_a], B[best_b]);
	value1 = A[best_a];
	value2 = B[best_b] * cudaSignum(dp);
}

//__global__ 
void cudaUpdateOrientation(int* phase, int num_phases, glm::dvec3* N, glm::dvec3* Q, Link* adj, int* adjOffset, int num_adj) {
//	int pi = blockIdx.x * blockDim.x + threadIdx.x;

	for (int pi = 0; pi < num_phases; ++pi) {
		if (pi >= num_phases)
			return;
		int i = phase[pi];
		glm::dvec3 n_i = N[i];
		double weight_sum = 0.0f;
		glm::dvec3 sum = Q[i];

		for (int l = adjOffset[i]; l < adjOffset[i + 1]; ++l) {
			Link link = adj[l];
			const int j = link.id;
			const double weight = link.weight;
			if (weight == 0)
				continue;
			glm::dvec3 n_j = N[j];
			glm::dvec3 q_j = Q[j];
			glm::dvec3 value1, value2;
			compat_orientation_extrinsic_4(sum, n_i, q_j, n_j, value1, value2);
			sum = value1 * weight_sum + value2 * weight;
			sum -= n_i*glm::dot(n_i, sum);
			weight_sum += weight;

			double norm = glm::length(sum);
			if (norm > 2.93873587705571876e-39f)
				sum /= norm;
		}

		if (weight_sum > 0) {
			Q[i] = sum;
		}
	}
}

//__global__
void cudaPropagateOrientationUpper(glm::dvec3* srcField, glm::ivec2* toUpper, glm::dvec3* N, glm::dvec3* destField, int num_orientation) {
//	int i = blockIdx.x * blockDim.x + threadIdx.x;
	for (int i = 0; i < num_orientation; ++i) {
		if (i >= num_orientation)
			return;
		for (int k = 0; k < 2; ++k) {
			int dest = toUpper[i][k];
			if (dest == -1)
				continue;
			glm::dvec3 q = srcField[i];
			glm::dvec3 n = N[dest];
			destField[dest] = q - n * glm::dot(n, q);
		}
	}
}

//__global__
void cudaPropagateOrientationLower(glm::ivec2* toUpper, glm::dvec3* Q, glm::dvec3* N, glm::dvec3* Q_next, glm::dvec3* N_next, int num_toUpper) {
//	int i = blockIdx.x * blockDim.x + threadIdx.x;
	for (int i = 0; i < num_toUpper; ++i) {
		if (i >= num_toUpper)
			return;
		glm::ivec2 upper = toUpper[i];
		glm::dvec3 q0 = Q[upper[0]];
		glm::dvec3 n0 = N[upper[0]];

		glm::dvec3 q, q1, n1, value1, value2;
		if (upper[1] != -1) {
			q1 = Q[upper[1]];
			n1 = N[upper[1]];
			compat_orientation_extrinsic_4(q0, n0, q1, n1, value1, value2);
			q = value1 + value2;
		}
		else {
			q = q0;
		}
		glm::dvec3 n = N_next[i];
		q -= glm::dot(n, q) * n;

		double len = q.x * q.x + q.y * q.y + q.z * q.z;
		if (len > 2.93873587705571876e-39f)
			q /= sqrt(len);
		Q_next[i] = q;
	}
}

void UpdateOrientation(int* phase, int num_phases, glm::dvec3* N, glm::dvec3* Q, Link* adj, int* adjOffset, int num_adj) {
//	cudaUpdateOrientation << <(num_phases + 255) / 256, 256 >> >(phase, num_phases, N, Q, adj, num_adj);
	cudaUpdateOrientation(phase, num_phases, N, Q, adj, adjOffset, num_adj);
}

void PropagateOrientationUpper(glm::dvec3* srcField, int num_orientation, glm::ivec2* toUpper, glm::dvec3* N, glm::dvec3* destField) {
//	cudaPropagateOrientationUpper << <(num_orientation + 255) / 256, 256 >> >(srcField, toUpper, N, destField, num_orientation);
	cudaPropagateOrientationUpper(srcField, toUpper, N, destField, num_orientation);
}

void PropagateOrientationLower(glm::ivec2* toUpper, glm::dvec3* Q, glm::dvec3* N, glm::dvec3* Q_next, glm::dvec3* N_next, int num_toUpper) {
//	cudaPropagateOrientationLower << <(num_toUpper + 255) / 256, 256 >> >(toUpper, Q, N, Q_next, num_toUpper);
	cudaPropagateOrientationLower(toUpper, Q, N, Q_next, N_next, num_toUpper);
}
