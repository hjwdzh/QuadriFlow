#ifndef FIELD_MATH_H_
#define FIELD_MATH_H_

#include <algorithm>
#include <Eigen/Core>
#include <Eigen/Dense>
using namespace Eigen;

inline float signum(float value) {
	return std::copysign((float)1, value);
}

/// Always-positive modulo function (assumes b > 0)
inline int modulo(int a, int b) {
	int r = a % b;
	return (r < 0) ? r + b : r;
}

inline std::pair<int, int>
compat_orientation_extrinsic_index_4(const Vector3f &q0, const Vector3f &n0,
const Vector3f &q1, const Vector3f &n1) {
	const Vector3f A[2] = { q0, n0.cross(q0) };
	const Vector3f B[2] = { q1, n1.cross(q1) };

	float best_score = -std::numeric_limits<float>::infinity();
	int best_a = 0, best_b = 0;

	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < 2; ++j) {
			float score = std::abs(A[i].dot(B[j]));
			if (score > best_score) {
				best_a = i;
				best_b = j;
				best_score = score;
			}
		}
	}

	if (A[best_a].dot(B[best_b]) < 0)
		best_b += 2;

	return std::make_pair(best_a, best_b);
}

inline std::pair<Vector3f, Vector3f>
compat_orientation_extrinsic_4(const Vector3f &q0, const Vector3f &n0,
const Vector3f &q1, const Vector3f &n1) {
	const Vector3f A[2] = { q0, n0.cross(q0) };
	const Vector3f B[2] = { q1, n1.cross(q1) };

	float best_score = -std::numeric_limits<float>::infinity();
	int best_a = 0, best_b = 0;

	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < 2; ++j) {
			float score = std::abs(A[i].dot(B[j]));
			if (score > best_score) {
				best_a = i;
				best_b = j;
				best_score = score;
			}
		}
	}

	const float dp = A[best_a].dot(B[best_b]);
	return std::make_pair(A[best_a], B[best_b] * signum(dp));
}


inline Vector3f middle_point(const Vector3f &p0, const Vector3f &n0, const Vector3f &p1, const Vector3f &n1) {
	/* How was this derived?
	*
	* Minimize \|x-p0\|^2 + \|x-p1\|^2, where
	* dot(n0, x) == dot(n0, p0)
	* dot(n1, x) == dot(n1, p1)
	*
	* -> Lagrange multipliers, set derivative = 0
	*  Use first 3 equalities to write x in terms of
	*  lambda_1 and lambda_2. Substitute that into the last
	*  two equations and solve for the lambdas. Finally,
	*  add a small epsilon term to avoid issues when n1=n2.
	*/
	float n0p0 = n0.dot(p0), n0p1 = n0.dot(p1),
		n1p0 = n1.dot(p0), n1p1 = n1.dot(p1),
		n0n1 = n0.dot(n1),
		denom = 1.0f / (1.0f - n0n1*n0n1 + 1e-4f),
		lambda_0 = 2.0f*(n0p1 - n0p0 - n0n1*(n1p0 - n1p1))*denom,
		lambda_1 = 2.0f*(n1p0 - n1p1 - n0n1*(n0p1 - n0p0))*denom;

	return 0.5f * (p0 + p1) - 0.25f * (n0 * lambda_0 + n1 * lambda_1);
}

inline Vector3f position_floor_4(const Vector3f &o, const Vector3f &q,
	const Vector3f &n, const Vector3f &p,
	float scale, float inv_scale) {
	Vector3f t = n.cross(q);
	Vector3f d = p - o;
	return o +
		q * std::floor(q.dot(d) * inv_scale) * scale +
		t * std::floor(t.dot(d) * inv_scale) * scale;
}

inline std::pair<Vector3f, Vector3f> compat_position_extrinsic_4(
	const Vector3f &p0, const Vector3f &n0, const Vector3f &q0, const Vector3f &o0,
	const Vector3f &p1, const Vector3f &n1, const Vector3f &q1, const Vector3f &o1,
	float scale, float inv_scale) {

	Vector3f t0 = n0.cross(q0), t1 = n1.cross(q1);
	Vector3f middle = middle_point(p0, n0, p1, n1);
	Vector3f o0p = position_floor_4(o0, q0, n0, middle, scale, inv_scale);
	Vector3f o1p = position_floor_4(o1, q1, n1, middle, scale, inv_scale);

	float best_cost = std::numeric_limits<float>::infinity();
	int best_i = -1, best_j = -1;

	for (int i = 0; i<4; ++i) {
		Vector3f o0t = o0p + (q0 * (i & 1) + t0 * ((i & 2) >> 1)) * scale;
		for (int j = 0; j<4; ++j) {
			Vector3f o1t = o1p + (q1 * (j & 1) + t1 * ((j & 2) >> 1)) * scale;
			float cost = (o0t - o1t).squaredNorm();

			if (cost < best_cost) {
				best_i = i;
				best_j = j;
				best_cost = cost;
			}
		}
	}

	return std::make_pair(
		o0p + (q0 * (best_i & 1) + t0 * ((best_i & 2) >> 1)) * scale,
		o1p + (q1 * (best_j & 1) + t1 * ((best_j & 2) >> 1)) * scale);
}

inline Vector3f position_round_4(const Vector3f &o, const Vector3f &q,
	const Vector3f &n, const Vector3f &p,
	float scale, float inv_scale) {
	Vector3f t = n.cross(q);
	Vector3f d = p - o;
	return o +
		q * std::round(q.dot(d) * inv_scale) * scale +
		t * std::round(t.dot(d) * inv_scale) * scale;
}

inline Vector2i position_floor_index_4(const Vector3f &o, const Vector3f &q,
	const Vector3f &n, const Vector3f &p,
	float /* unused */, float inv_scale) {
	Vector3f t = n.cross(q);
	Vector3f d = p - o;
	return Vector2i(
		(int)std::floor(q.dot(d) * inv_scale),
		(int)std::floor(t.dot(d) * inv_scale));
}

inline std::pair<Vector2i, Vector2i> compat_position_extrinsic_index_4(
	const Vector3f &p0, const Vector3f &n0, const Vector3f &q0, const Vector3f &o0,
	const Vector3f &p1, const Vector3f &n1, const Vector3f &q1, const Vector3f &o1,
	float scale, float inv_scale, float* error) {
	Vector3f t0 = n0.cross(q0), t1 = n1.cross(q1);
	Vector3f middle = middle_point(p0, n0, p1, n1);
	Vector2i o0p = position_floor_index_4(o0, q0, n0, middle, scale, inv_scale);
	Vector2i o1p = position_floor_index_4(o1, q1, n1, middle, scale, inv_scale);

	float best_cost = std::numeric_limits<float>::infinity();
	int best_i = -1, best_j = -1;

	for (int i = 0; i<4; ++i) {
		Vector3f o0t = o0 + (q0 * ((i & 1) + o0p[0]) + t0 * (((i & 2) >> 1) + o0p[1])) * scale;
		for (int j = 0; j<4; ++j) {
			Vector3f o1t = o1 + (q1 * ((j & 1) + o1p[0]) + t1 * (((j & 2) >> 1) + o1p[1])) * scale;
			float cost = (o0t - o1t).squaredNorm();

			if (cost < best_cost) {
				best_i = i;
				best_j = j;
				best_cost = cost;
			}
		}
	}
	if (error)
		*error = best_cost;

	return std::make_pair(
		Vector2i((best_i & 1) + o0p[0], ((best_i & 2) >> 1) + o0p[1]),
		Vector2i((best_j & 1) + o1p[0], ((best_j & 2) >> 1) + o1p[1]));
}

inline void coordinate_system(const Vector3f &a, Vector3f &b, Vector3f &c) {
	if (std::abs(a.x()) > std::abs(a.y())) {
		float invLen = 1.0f / std::sqrt(a.x() * a.x() + a.z() * a.z());
		c = Vector3f(a.z() * invLen, 0.0f, -a.x() * invLen);
	}
	else {
		float invLen = 1.0f / std::sqrt(a.y() * a.y() + a.z() * a.z());
		c = Vector3f(0.0f, a.z() * invLen, -a.y() * invLen);
	}
	b = c.cross(a);
}

inline Vector3f rotate_vector_into_plane(Vector3f q, const Vector3f &source_normal, const Vector3f &target_normal) {
	const float cosTheta = source_normal.dot(target_normal);
	if (cosTheta < 0.9999f) {
		Vector3f axis = source_normal.cross(target_normal);
		q = q * cosTheta + axis.cross(q) +
			axis * (axis.dot(q) * (1.0f - cosTheta) / axis.dot(axis));
	}
	return q;
}

inline Vector3f Travel(Vector3f p, const Vector3f& dir, float& len, int& f, VectorXi& E2E, MatrixXf& V, MatrixXi& F, MatrixXf& NF) {
	Vector3f N = NF.col(f);
	Vector3f pt = (dir - dir.dot(N) * N).normalized();
	int prev_id = -1;
	int count = 0;
	while (len > 0) {
		count += 1;
		Vector3f N = NF.col(f);
		int edge_id = f * 3;
		double max_len = 1e30;
		bool found = false;
		int next_id, next_f;
		for (int fid = 0; fid < 3; ++fid) {
			if (fid + edge_id == prev_id)
				continue;
			Vector3f dir1 = V.col(F(fid, f)) - p;
			Vector3f dir2 = V.col(F((fid + 1) % 3, f)) - p;
			Vector3f q = (dir2 - dir1).normalized();
			Vector3f h = dir1 - dir1.dot(q) * q;
			float h_sum = h.norm();
			h /= h_sum;
			float unit_h = pt.dot(h);
			float t = h_sum / unit_h;

			Vector3f qq = pt * t;
			for (int j = 0; j < 3; ++j) {
				printf("%f ", (qq[j] - dir1[j]) / (dir2[j] - dir1[j]));
			}
			printf("\n");

			if (t >= 0 && t < max_len) {
				printf("#############\n");
				max_len = t;
				next_id = E2E[edge_id + fid];
				next_f = next_id;
				if (next_f != -1)
					next_f /= 3;
				found = true;
			}
		}
		if (!found) {
			printf("error...\n");
			system("pause");
		}
		printf("%f %f %d\n", len, max_len, f);
		if (max_len >= len) {
			p = p + len * pt;
			len = 0;
			return p;
		}
		p = p + max_len * pt;
		len -= max_len;
		if (next_f == -1)
			return p;
		pt = rotate_vector_into_plane(pt, NF.col(f), NF.col(next_f));
		f = next_f;
		prev_id = next_id;

		int pid = prev_id % 3;
		Vector3f v1 = V.col(F(pid, f));
		Vector3f v2 = V.col(F((pid + 1) % 3, f));
		printf("%f %f %f\n", (p.x() - v1.x()) / (v2.x() - v1.x()),
			(p.y() - v1.y()) / (v2.y() - v1.y()),
			(p.z() - v1.z()) / (v2.z() - v1.z()));
		pid = pid;
	}
	return p;
}

#endif