#include "Optimizer.h"
#include "field_math.h"
#include <fstream>

Optimizer::Optimizer()
{}

void Optimizer::optimize_orientations(Hierarchy &mRes)
{
	int levelIterations = 6;
	for (int level = mRes.mAdj.size() - 1; level >= 0; --level) {
		AdjacentMatrix &adj = mRes.mAdj[level];
		const MatrixXf &N = mRes.mN[level];
		MatrixXf &Q = mRes.mQ[level];

		for (int iter = 0; iter < levelIterations; ++iter) {
			for (int i = 0; i < N.cols(); ++i) {
				const Vector3f n_i = N.col(i);
				float weight_sum = 0.0f;
				Vector3f sum = Q.col(i);
				for (auto& link : adj[i]) {
					const int j = link.id;
					const float weight = link.weight;
					if (weight == 0)
						continue;
					const Vector3f n_j = N.col(j);
					Vector3f q_j = Q.col(j);
					std::pair<Vector3f, Vector3f> value = compat_orientation_extrinsic_4(sum, n_i, q_j, n_j);
					sum = value.first * weight_sum + value.second * weight;
					sum -= n_i*n_i.dot(sum);
					weight_sum += weight;

					float norm = sum.norm();
					if (norm > RCPOVERFLOW)
						sum /= norm;
				}

				if (weight_sum > 0)
					Q.col(i) = sum;
			}
		}

		if (level > 0) {
			const MatrixXf &srcField = mRes.mQ[level];
			const MatrixXi &toUpper = mRes.mToUpper[level - 1];
			MatrixXf &destField = mRes.mQ[level - 1];
			const MatrixXf &N = mRes.mN[level - 1];
			for (int i = 0; i < srcField.cols(); ++i) {
				for (int k = 0; k < 2; ++k) {
					int dest = toUpper(k, i);
					if (dest == -1)
						continue;
					Vector3f q = srcField.col(i), n = N.col(dest);
					destField.col(dest) = q - n * n.dot(q);
				}
			}
		}
	}
	for (int l = 0; l< mRes.mN.size() - 1; ++l)  {
		const MatrixXf &N = mRes.mN[l];
		const MatrixXf &N_next = mRes.mN[l + 1];
		const MatrixXf &Q = mRes.mQ[l];
		MatrixXf &Q_next = mRes.mQ[l + 1];
		auto& toUpper = mRes.mToUpper[l];
		for (int i = 0; i < toUpper.cols(); ++i) {
			Vector2i upper = toUpper.col(i);
			Vector3f q0 = Q.col(upper[0]);
			Vector3f n0 = N.col(upper[0]);
			Vector3f q;

			if (upper[1] != -1) {
				Vector3f q1 = Q.col(upper[1]);
				Vector3f n1 = N.col(upper[1]);
				auto result = compat_orientation_extrinsic_4(q0, n0, q1, n1);
				q = result.first + result.second;
			}
			else {
				q = q0;
			}
			Vector3f n = N_next.col(i);
			q -= n.dot(q) * n;
			if (q.squaredNorm() > RCPOVERFLOW)
				q.normalize();

			Q_next.col(i) = q;
		}
	}
}

void Optimizer::optimize_scale(Hierarchy &mRes)
{
	int levelIterations = 6;
	for (int level = mRes.mAdj.size() - 1; level >= 0; --level) {
		if (level != 0)
			continue;
		AdjacentMatrix &adj = mRes.mAdj[level];
		const MatrixXf &N = mRes.mN[level];
		MatrixXf &Q = mRes.mQ[level];
		MatrixXf &V = mRes.mV[level];
		MatrixXf &S = mRes.mS[level];
		MatrixXf dscale(2, mRes.mV[level].cols());
		for (int iter = -1; iter < levelIterations; ++iter) {
			for (int i = 0; i < N.cols(); ++i) {
				const Vector3f n_i = N.col(i);
				const Vector3f p_i = V.col(i);
				float weight_sum[2] = { 0.0f, 0.0f };
				Vector2f dsum(0.0f, 0.0f);
				Vector3f q_i = Q.col(i);
				std::vector<std::pair<float, std::pair<float, float> > > buffer[2];
				for (auto& link : adj[i]) {
					const int j = link.id;
					const float weight = link.weight;
					if (weight == 0)
						continue;
					const Vector3f n_j = N.col(j);
					Vector3f q_j = Q.col(j);
					q_j = rotate_vector_into_plane(q_j, n_j, n_i);

					const Vector3f p_j = V.col(j);
					Vector3f A[2] = { q_i, n_i.cross(q_i) };
					Vector3f B[2] = { q_j, n_i.cross(q_j) };

					float best_score = -std::numeric_limits<float>::infinity();
					int best_a = 0, best_b = 0;
					for (int l = 0; l < 2; ++l) {
						for (int k = 0; k < 2; ++k) {
							float score = std::abs(A[l].dot(B[k]));
							if (score > best_score) {
								best_a = l;
								best_b = k;
								best_score = score;
							}
						}
					}
					if (best_a == 1) {
						best_a = 0;
						best_b = 1 - best_b;
					}
					float sign[2] = { 1, 1 };
					if (A[0].dot(B[best_b]) < 0)
						sign[0] = -1;
					if (A[1].dot(B[1 - best_b]) < 0)
						sign[1] = -1;
					if (iter == -1) {
						for (int k = 0; k < 2; ++k) {
							int b_ind = 1 - (best_b + k) % 2;
							double dis = (p_j - p_i).dot(A[k]);
							double diff = B[b_ind].dot(A[k]) * sign[k];
							if (fabs(dis) <= RCPOVERFLOW)
								continue;
							buffer[k].push_back(std::make_pair(fabs(dis), std::make_pair(diff / dis, link.weight)));
						}
					}
					else {
						for (int k = 0; k < 2; ++k) {
							int b_ind = 1 - (best_b + k) % 2;
							double dis = (p_i - p_j).dot(B[b_ind]);
							if (fabs(dis) <= RCPOVERFLOW)
								continue;
							double new_s = S(b_ind, j) + dscale(k, j) * dis * sign[k];
							dsum(k) = new_s * link.weight + dsum(k) * weight_sum[k];
							weight_sum[k] += link.weight;
							if (weight_sum[k] > RCPOVERFLOW)
								dsum(k) /= weight_sum[k];
						}
					}
				}
				if (iter == -1) {
					for (int k = 0; k < 2; ++k) {
						std::sort(buffer[k].rbegin(), buffer[k].rend());
						float sum_1 = 0.0f, sum_2 = 0.0f;
						for (int j = 0; j < buffer[k].size(); ++j) {
							if (buffer[k][j].first < buffer[k][0].first * 0.3)
								break;
							sum_1 += buffer[k][j].second.first * buffer[k][j].second.second;
							sum_2 += buffer[k][j].second.second;
						}
						dscale(k, i) = sum_1 / sum_2;
					}
				}
				else {
					S.col(i) = dscale.col(i);
				}
			}
		}

		if (level > 0) {
			const MatrixXf &srcField = mRes.mS[level];
			const MatrixXi &toUpper = mRes.mToUpper[level - 1];
			MatrixXf &destField = mRes.mS[level - 1];
			float scale_sum = 0.0f;
			for (int i = 0; i < srcField.cols(); ++i) {
				for (int k = 0; k < 2; ++k) {
					scale_sum += srcField(k, i);
				}
			}
			scale_sum = (2.0f * srcField.cols()) / scale_sum;
			for (int i = 0; i < srcField.cols(); ++i) {
				for (int k = 0; k < 2; ++k) {
					int dest = toUpper(k, i);
					if (dest == -1)
						continue;
					Vector2f q = srcField.col(i);
					destField.col(dest) = q * scale_sum;
				}
			}
		}
		else {
			MatrixXf &srcField = mRes.mS[level];
			float maxf[2] = { -1e30, -1e30 };
			float minf[2] = { 1e30, 1e30 };
			for (int i = 0; i < srcField.cols(); ++i) {
				for (int k = 0; k < 2; ++k) {
					maxf[k] = fmax(maxf[k], srcField(k, i));
					minf[k] = fmin(minf[k], srcField(k, i));
				}
			}
			printf("%f %f %f %f\n", minf[0], maxf[0], minf[1], maxf[1]);
			for (int i = 0; i < srcField.cols(); ++i) {
				for (int k = 0; k < 2; ++k) {
					float s = (srcField(k, i) + 5) / 10.0f;//(maxf[k] - minf[k]);
//					if (k == 0)
//					printf("%f %f\n", srcField(k, i), s);
					if (s < 0)
						s = 0;
					if (s > 1)
						s = 1;
					srcField(k, i) = s;
				}
			}
		}
	}
}

void Optimizer::optimize_positions(Hierarchy &mRes)
{
	int levelIterations = 6;

	for (int level = mRes.mAdj.size() - 1; level >= 0; --level) {
		AdjacentMatrix &adj = mRes.mAdj[level];
		const MatrixXf &N = mRes.mN[level];
		MatrixXf &Q = mRes.mQ[level];

		for (int iter = 0; iter < levelIterations; ++iter) {
			AdjacentMatrix &adj = mRes.mAdj[level];
			const MatrixXf &N = mRes.mN[level], &Q = mRes.mQ[level], &V = mRes.mV[level];
			const float scale = mRes.mScale, inv_scale = 1.0f / scale;
			MatrixXf &O = mRes.mO[level];

			for (int i = 0; i < N.cols(); ++i) {
				const Vector3f n_i = N.col(i), v_i = V.col(i);
				Vector3f q_i = Q.col(i);

				Vector3f sum = O.col(i);
				float weight_sum = 0.0f;

				q_i.normalize();
				for (auto& link : adj[i]) {
					const int j = link.id;
					const float weight = link.weight;
					if (weight == 0)
						continue;

					const Vector3f n_j = N.col(j), v_j = V.col(j);
					Vector3f q_j = Q.col(j), o_j = O.col(j);

					q_j.normalize();

					std::pair<Vector3f, Vector3f> value = compat_position_extrinsic_4(
						v_i, n_i, q_i, sum, v_j, n_j, q_j, o_j, scale, inv_scale);

					sum = value.first*weight_sum + value.second*weight;
					weight_sum += weight;
					if (weight_sum > RCPOVERFLOW)
						sum /= weight_sum;
					sum -= n_i.dot(sum - v_i)*n_i;
				}

				if (weight_sum > 0) {
					O.col(i) = position_round_4(sum, q_i, n_i, v_i, scale, inv_scale);
				}
			}
		}

		if (level > 0) {
			const MatrixXf &srcField = mRes.mO[level];
			const MatrixXi &toUpper = mRes.mToUpper[level - 1];
			MatrixXf &destField = mRes.mO[level - 1];
			const MatrixXf &N = mRes.mN[level - 1];
			const MatrixXf &V = mRes.mV[level - 1];
			for (int i = 0; i < srcField.cols(); ++i) {
				for (int k = 0; k < 2; ++k) {
					int dest = toUpper(k, i);
					if (dest == -1)
						continue;
					Vector3f o = srcField.col(i), n = N.col(dest), v = V.col(dest);
					o -= n * n.dot(o - v);
					destField.col(dest) = o;
				}
			}
		}
	}
}

