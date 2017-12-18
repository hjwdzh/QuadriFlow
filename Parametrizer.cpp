#include "Parametrizer.h"

#include "loader.h"
#include "MergeVertex.h"
#include "subdivide.h"
#include "dedge.h"
#include "AdjacentMatrix.h"
#include "field_math.h"

#include <fstream>

#define LOG_OUTPUT

inline float fast_acos(float x) {
	float negate = float(x < 0.0f);
	x = std::abs(x);
	float ret = -0.0187293f;
	ret *= x; ret = ret + 0.0742610f;
	ret *= x; ret = ret - 0.2121144f;
	ret *= x; ret = ret + 1.5707288f;
	ret = ret * std::sqrt(1.0f - x);
	ret = ret - 2.0f * negate * ret;
	return negate * (float)M_PI + ret;
}

void Parametrizer::Load(const char* filename)
{
	load(filename, V, F);
#ifdef LOG_OUTPUT
	printf("vertices size: %d\n", V.cols());
	printf("faces size: %d\n", F.cols());
#endif

	merge_close(V, F, 1e-6);

}

void Parametrizer::Initialize()
{
	ComputeMeshStatus();

	num_vertices = V.cols() / 16;
	num_faces = num_vertices;
	scale = sqrt(surface_area / num_faces);

	printf("Compute Direct Graph\n");
	compute_direct_graph(V, F, V2E, E2E, boundary, nonManifold);
	printf("Compute Direct Graph finish\n");

	float target_len = std::min(scale / 2, average_edge_length * 2);
	if (target_len < max_edge_length) {
		subdivide(F, V, V2E, E2E, boundary, nonManifold, target_len);
		compute_direct_graph(V, F, V2E, E2E, boundary, nonManifold);
	}
	printf("Adjacency Matrix\n");
	generate_adjacency_matrix_uniform(F, V2E, E2E, nonManifold, adj);
	printf("Adjacency Matrix finish\n");
	ComputeSmoothNormal();
	ComputeVertexArea();

	hierarchy.mA[0] = std::move(A);
	hierarchy.mAdj[0] = std::move(adj);
	hierarchy.mN[0] = std::move(N);
	hierarchy.mV[0] = std::move(V);
	hierarchy.mE2E = std::move(E2E);
	hierarchy.mF = std::move(F);
	hierarchy.Initialize(scale);
}

void Parametrizer::ComputeMeshStatus()
{
	surface_area = 0;
	average_edge_length = 0;
	max_edge_length = 0;
	for (int f = 0; f < F.cols(); ++f) {
		Vector3f v[3] = { V.col(F(0, f)), V.col(F(1, f)), V.col(F(2, f)) };
		float area = 0.5f * (v[1] - v[0]).cross(v[2] - v[0]).norm();
		surface_area += area;
		for (int i = 0; i < 3; ++i) {
			float len = (v[(i + 1) % 3] - v[i]).norm();
			average_edge_length += len;
			if (len > max_edge_length)
				max_edge_length = len;
		}
	}
	average_edge_length /= (F.cols() * 3);
}

void Parametrizer::ComputeSmoothNormal()
{
	/* Compute face normals */
	MatrixXf Nf(3, F.cols());
	for (int f = 0; f < F.cols(); ++f) {
		Vector3f v0 = V.col(F(0, f)),
			v1 = V.col(F(1, f)),
			v2 = V.col(F(2, f)),
			n = (v1 - v0).cross(v2 - v0);
		float norm = n.norm();
		if (norm < RCPOVERFLOW) {
			n = Vector3f::UnitX();
		}
		else {
			n /= norm;
		}
		Nf.col(f) = n;
	}

	N.resize(3, V.cols());
	for (int i = 0; i < V2E.rows(); ++i) {
		int edge = V2E[i];
		if (nonManifold[i] || edge == -1) {
			N.col(i) = Vector3f::UnitX();
			continue;
		}

		int stop = edge;
		Vector3f normal = Vector3f::Zero();
		do {
			int idx = edge % 3;

			Vector3f d0 = V.col(F((idx + 1) % 3, edge / 3)) - V.col(i);
			Vector3f d1 = V.col(F((idx + 2) % 3, edge / 3)) - V.col(i);
			float angle = fast_acos(d0.dot(d1) / std::sqrt(d0.squaredNorm() * d1.squaredNorm()));

			/* "Computing Vertex Normals from Polygonal Facets"
			by Grit Thuermer and Charles A. Wuethrich, JGT 1998, Vol 3 */
			if (std::isfinite(angle))
				normal += Nf.col(edge / 3) * angle;

			int opp = E2E[edge];
			if (opp == -1)
				break;

			edge = dedge_next_3(opp);
		} while (edge != stop);
		float norm = normal.norm();
		N.col(i) = norm > RCPOVERFLOW ? Vector3f(normal / norm)
			: Vector3f::UnitX();
	}

}


void Parametrizer::ComputeVertexArea() {
	A.resize(V.cols());
	A.setZero();

	for (int i = 0; i < V2E.size(); ++i) {
		int edge = V2E[i], stop = edge;
		if (nonManifold[i] || edge == -1)
			continue;
		float vertex_area = 0;
		do {
			int ep = dedge_prev_3(edge), en = dedge_next_3(edge);

			Vector3f v = V.col(F(edge % 3, edge / 3));
			Vector3f vn = V.col(F(en % 3, en / 3));
			Vector3f vp = V.col(F(ep % 3, ep / 3));

			Vector3f face_center = (v + vp + vn) * (1.0f / 3.0f);
			Vector3f prev = (v + vp) * 0.5f;
			Vector3f next = (v + vn) * 0.5f;

			vertex_area +=
				0.5f * ((v - prev).cross(v - face_center).norm() +
				(v - next).cross(v - face_center).norm());

			int opp = E2E[edge];
			if (opp == -1)
				break;
			edge = dedge_next_3(opp);
		} while (edge != stop);

		A[i] = vertex_area;
	}


}

#include <set>
#include <map>
#include "dset.h"

void extract_graph(const Hierarchy &mRes,
	std::vector<std::vector<TaggedLink> > &adj_new,
	MatrixXf &O_new, MatrixXf &N_new)
{
	int numV = mRes.mV[0].cols();
	DisjointSets dset(numV);
	adj_new.clear();
	adj_new.resize(numV);

	typedef std::pair<int, int> Edge;
	typedef std::pair<Edge, float> WeightedEdge;
	std::vector<WeightedEdge> collapse_edge_vec;
	collapse_edge_vec.reserve((uint32_t)(numV*2.5f));

	const MatrixXf &Q = mRes.mQ[0], &O = mRes.mO[0], &N = mRes.mN[0], &V = mRes.mV[0];
	const float scale = mRes.mScale, inv_scale = 1.0f / scale;

	// clasify
	for (int i = 0; i < numV; ++i) {
		for (auto& link : mRes.mAdj[0][i]) {
			int j = link.id;

			if (j < i)
				continue;

			std::pair<Vector3f, Vector3f> Q_rot = compat_orientation_extrinsic_4(
				Q.col(i), N.col(i), Q.col(j), N.col(j));

			float error = 0;
			const auto Vi = V.col(i), Ni = N.col(i), Oi = O.col(i);
			const auto Vj = V.col(j), Nj = N.col(j), Oj = O.col(j);
			std::pair<Vector2i, Vector2i> shift = compat_position_extrinsic_index_4(
				V.col(i), N.col(i), Q_rot.first, O.col(i),
				V.col(j), N.col(j), Q_rot.second, O.col(j),
				scale, inv_scale, &error);

			Vector2i absDiff = (shift.first - shift.second).cwiseAbs();

			if (absDiff.maxCoeff() > 1 || (absDiff == Vector2i(1, 1)))
				continue; /* Ignore longer-distance links and diagonal lines for quads */
			bool collapse = absDiff.sum() == 0;

			if (collapse) {
				collapse_edge_vec.push_back(std::make_pair(std::make_pair(i, j), error));
			}
			else {
				adj_new[i].push_back(j);
				adj_new[j].push_back(i);
			}
		}
	}
	// collapse neighbors
	struct WeightedEdgeComparator {
		bool operator()(const WeightedEdge& e1, const WeightedEdge& e2) const { return e1.second < e2.second; }
	};
	std::stable_sort(collapse_edge_vec.begin(), collapse_edge_vec.end());
	std::vector<int> nCollapses(numV, 0);
	int nItem = 0;
	for (int i = 0; i < collapse_edge_vec.size(); ++i) {
		const WeightedEdge &we = collapse_edge_vec[nItem++];
		Edge edge = we.first;
		int t1 = dset.find(edge.first);
		int t2 = dset.find(edge.second);
		if (t1 == t2)
			continue;
		int target = dset.unite(t1, t2);
		std::set<int> temp;
		for (auto& t : adj_new[t1]) {
			temp.insert(t.id);
		}
		for (auto& t : adj_new[t2]) {
			temp.insert(t.id);
		}
		adj_new[t1].clear();
		adj_new[t2].clear();
		nCollapses[target] = nCollapses[t1] + nCollapses[t2] + 1;
		adj_new[target].insert(adj_new[target].end(), temp.begin(), temp.end());
	}
	// compress the array
	int nVertices = 0;
	std::map<int, int> vertex_map;
	float avg_collapses = 0;

	for (uint32_t i = 0; i<adj_new.size(); ++i) {
		if (adj_new[i].empty())
			continue;
		if (i != nVertices) {
			adj_new[nVertices].swap(adj_new[i]);
			std::swap(nCollapses[nVertices], nCollapses[i]);
		}
		avg_collapses += nCollapses[nVertices];
		vertex_map[i] = nVertices++;
	}
	adj_new.resize(nVertices);
	adj_new.shrink_to_fit();
	avg_collapses /= nVertices;

	for (int i = 0; i < adj_new.size(); ++i) {
		std::set<int> temp;
		for (auto& k : adj_new[i])
			temp.insert(vertex_map[dset.find(k.id)]);
		std::vector<TaggedLink> new_vec;
		new_vec.reserve(temp.size());
		for (auto& j : temp)
			new_vec.push_back(TaggedLink(j));
		adj_new[i] = std::move(new_vec);
	}

	// Remove Spurious Vertices
	int removed = 0;
	for (int i = 0; i<adj_new.size(); ++i) {
		if (nCollapses[i] > avg_collapses / 10)
			continue;

		for (auto neighbor : adj_new[i]) {
			auto &a = adj_new[neighbor.id];
			a.erase(std::remove_if(a.begin(), a.end(), [&](const TaggedLink &v) { return v.id == i; }), a.end());
		}

		adj_new[i].clear();
		++removed;
	}

	// Compute O_new and N_new
	O_new.resize(3, nVertices);
	N_new.resize(3, nVertices);
	O_new.setZero();
	N_new.setZero();

	Eigen::VectorXf cluster_weight(nVertices);
	cluster_weight.setZero();

	for (int i = 0; i < numV; ++i) {
		auto it = vertex_map.find(dset.find(i));
		if (it == vertex_map.end())
			continue;
		uint32_t j = it->second;

		float weight = std::exp(-(O.col(i) - V.col(i)).squaredNorm() * inv_scale * inv_scale * 9);
		for (uint32_t k = 0; k < 3; ++k) {
			O_new.coeffRef(k, j) += O(k, i) * weight;
			N_new.coeffRef(k, j) += N(k, i) * weight;
		}
		cluster_weight[j] += weight;
	}
	for (uint32_t i = 0; i<nVertices; ++i) {
		if (cluster_weight[i] == 0) {
			printf("waring %d %d\n", i, adj_new[i].size());
			system("pause");
			continue;
		}
		O_new.col(i) /= cluster_weight[i];
		N_new.col(i).normalize();
	}

	// remove unnecessary edges...
	// sort edges
	for (int i = 0; i < O_new.cols(); ++i) {
		Vector3f s, t, p = O_new.col(i);
		coordinate_system(N_new.col(i), s, t);

		std::sort(adj_new[i].begin(), adj_new[i].end(),
			[&](const TaggedLink &j0, const TaggedLink &j1) {
			Vector3f v0 = O_new.col(j0.id) - p, v1 = O_new.col(j1.id) - p;
			return std::atan2(t.dot(v0), s.dot(v0)) > std::atan2(t.dot(v1), s.dot(v1));
		});
	}

	// remove unnecessary edges
	bool changed;
	uint32_t nRemoved = 0, nSnapped = 0;
	do {
		changed = false;

		bool changed_inner;
		do {
			changed_inner = false;
			float thresh = 0.3f * scale;

			std::vector<std::tuple<float, uint32_t, uint32_t, uint32_t>> candidates;
			for (uint32_t i_id = 0; i_id<adj_new.size(); ++i_id) {
				auto const &adj_i = adj_new[i_id];
				const Vector3f p_i = O_new.col(i_id);
				for (uint32_t j = 0; j<adj_i.size(); ++j) {
					uint32_t j_id = adj_i[j].id;
					const Vector3f p_j = O_new.col(j_id);
					auto const &adj_j = adj_new[j_id];

					for (uint32_t k = 0; k<adj_j.size(); ++k) {
						uint32_t k_id = adj_j[k].id;
						if (k_id == i_id)
							continue;
						const Vector3f p_k = O_new.col(k_id);
						float a = (p_j - p_k).norm(), b = (p_i - p_j).norm(), c = (p_i - p_k).norm();
						if (a > std::max(b, c)) {
							float s = 0.5f * (a + b + c);
							float height = 2 * std::sqrt(s*(s - a)*(s - b)*(s - c)) / a;
							if (height < thresh)
								candidates.push_back(std::make_tuple(height, i_id, j_id, k_id));
						}
					}
				}
			}

			std::sort(candidates.begin(), candidates.end(), [&](
				const decltype(candidates)::value_type &v0,
				const decltype(candidates)::value_type &v1)
			{ return std::get<0>(v0) < std::get<0>(v1); });

			for (auto t : candidates) {
				uint32_t i = std::get<1>(t), j = std::get<2>(t), k = std::get<3>(t);
				bool edge1 = std::find_if(adj_new[i].begin(), adj_new[i].end(),
					[j](const TaggedLink &l) { return l.id == j; }) != adj_new[i].end();
				bool edge2 = std::find_if(adj_new[j].begin(), adj_new[j].end(),
					[k](const TaggedLink &l) { return l.id == k; }) != adj_new[j].end();
				bool edge3 = std::find_if(adj_new[k].begin(), adj_new[k].end(),
					[i](const TaggedLink &l) { return l.id == i; }) != adj_new[k].end();

				if (!edge1 || !edge2)
					continue;

				const Vector3f p_i = O_new.col(i), p_j = O_new.col(j), p_k = O_new.col(k);
				float a = (p_j - p_k).norm(), b = (p_i - p_j).norm(), c = (p_i - p_k).norm();
				float s = 0.5f * (a + b + c);
				float height = 2 * std::sqrt(s*(s - a)*(s - b)*(s - c)) / a;
				if (height != std::get<0>(t))
					continue;
				if ((p_i - p_j).norm() < thresh || (p_i - p_k).norm() < thresh) {
					uint32_t merge_id = (p_i - p_j).norm() < thresh ? j : k;
					O_new.col(i) = (O_new.col(i) + O_new.col(merge_id)) * 0.5f;
					N_new.col(i) = (N_new.col(i) + N_new.col(merge_id)) * 0.5f;
					std::set<uint32_t> adj_updated;
					for (auto const &n : adj_new[merge_id]) {
						if (n.id == i)
							continue;
						adj_updated.insert(n.id);
						for (auto &n2 : adj_new[n.id]) {
							if (n2.id == merge_id)
								n2.id = i;
						}
					}
					for (auto &n : adj_new[i])
						adj_updated.insert(n.id);
					adj_updated.erase(i);
					adj_updated.erase(merge_id);
					adj_new[merge_id].clear();
					adj_new[i].clear();
					for (uint32_t l : adj_updated)
						adj_new[i].push_back(l);
				}
				else {
					Vector3f n_k = N_new.col(k), n_j = N_new.col(j);
					//Vector3f dp = p_k - p_j, dn = n_k - n_j;
					//Float t = dp.dot(p_i-p_j) / dp.dot(dp);
					//O_new.col(i) = p_j + t*dp;
					//N_new.col(i) = (n_j + t*dn).normalized();
					O_new.col(i) = (p_j + p_k) * 0.5f;
					N_new.col(i) = (n_j + n_k).normalized();

					adj_new[j].erase(std::remove_if(adj_new[j].begin(), adj_new[j].end(),
						[k](const TaggedLink &l) { return l.id == k; }), adj_new[j].end());
					adj_new[k].erase(std::remove_if(adj_new[k].begin(), adj_new[k].end(),
						[j](const TaggedLink &l) { return l.id == j; }), adj_new[k].end());

					if (!edge3) {
						adj_new[i].push_back(k);
						adj_new[k].push_back(i);
					}
				}

				changed = true;
				changed_inner = true;
				++nSnapped;
			}
		} while (changed_inner);

		std::vector<std::pair<float, Edge>> candidates;
		for (size_t i = 0; i<adj_new.size(); ++i) {
			auto const &adj_i = adj_new[i];
			const Vector3f p_i = O_new.col(i);
			for (uint32_t j = 0; j<adj_i.size(); ++j) {
				uint32_t j_id = adj_i[j].id;
				const Vector3f p_j = O_new.col(j_id);

				uint32_t nTris = 0;
				float length = 0.0f;
				for (uint32_t k = 0; k<adj_i.size(); ++k) {
					uint32_t k_id = adj_i[k].id;
					if (k_id == j_id)
						continue;
					const Vector3f p_k = O_new.col(k_id);
					if (std::find_if(adj_new[j_id].begin(), adj_new[j_id].end(),
						[k_id](const TaggedLink &l) { return l.id == k_id; }) == adj_new[j_id].end())
						continue;
					nTris++;
					length += (p_k - p_i).norm() + (p_k - p_j).norm();
				}

				if (nTris == 2) {
					float exp_diag = length / 4 * std::sqrt(2.f);
					float diag = (p_i - p_j).norm();
					float score = std::abs((diag - exp_diag) / std::min(diag, exp_diag));
					candidates.push_back(std::make_pair(std::abs(score), std::make_pair(i, j_id)));
				}
			}
			std::sort(candidates.begin(), candidates.end(), [&](
				const std::pair<float, Edge> &v0,
				const std::pair<float, Edge> &v1) { return v0.first < v1.first; });

			for (auto c : candidates) {
				uint32_t i_id = c.second.first, j_id = c.second.second;
				auto const &adj_i = adj_new[i_id];
				uint32_t nTris = 0;
				for (uint32_t k = 0; k<adj_i.size(); ++k) {
					uint32_t k_id = adj_i[k].id;
					if (std::find_if(adj_new[j_id].begin(), adj_new[j_id].end(),
						[k_id](const TaggedLink &l) { return l.id == k_id; }) == adj_new[j_id].end())
						continue;
					nTris++;
				}
				if (nTris == 2) {
					adj_new[i_id].erase(std::remove_if(adj_new[i_id].begin(), adj_new[i_id].end(),
						[j_id](const TaggedLink &l) { return l.id == j_id; }), adj_new[i_id].end());
					adj_new[j_id].erase(std::remove_if(adj_new[j_id].begin(), adj_new[j_id].end(),
						[i_id](const TaggedLink &l) { return l.id == i_id; }), adj_new[j_id].end());
					changed = true;
					++nRemoved;
				}
			}
		}
	} while (changed);

}


void Parametrizer::ExtractMesh() {
	extract_graph(hierarchy, adj_extracted,
		mV_extracted, mN_extracted);
	/*
	Vector3f red = Vector3f::UnitX();

	int smooth_iterations = (int)(mSmoothSlider->value() * 10);
	extract_faces(adj_extracted, mV_extracted, mN_extracted, mNf_extracted,
		mF_extracted, posy, mRes.scale(), creaseOut, true,
		mPureQuadBox->checked(), mBVH, smooth_iterations);

	cout << "Extraction is done. (total time: " << timeString(timer.value()) << ")" << endl;

	int fmult = posy == 3 ? 1 : 2;

	MatrixXu F_gpu(3, mF_extracted.cols()*fmult);
	MatrixXf N_gpu(3, mF_extracted.cols()*posy);
	MatrixXf O_gpu(3, mF_extracted.cols()*posy);
	MatrixXf outputMeshWireframe(6, mF_extracted.cols() * posy * 2);

	for (uint32_t i = 0; i<(uint32_t)mF_extracted.cols(); ++i) {
		int base = posy*i;
		F_gpu.col(fmult*i) = Vector3u(base + 0, base + 1, base + 2);
		if (posy == 4)
			F_gpu.col(fmult*i + 1) = Vector3u(base + 2, base + 3, base + 0);
		bool irregular = posy == 4 && mF_extracted(2, i) == mF_extracted(3, i);
		for (int j = 0; j<posy; ++j) {
			uint32_t k = mF_extracted(j, i), kn = mF_extracted((j + 1) % posy, i);

			Vector3f col = red;
			if (irregular && j >= 1)
				col = Vector3f::Zero();
			outputMeshWireframe.col(i * 2 * posy + j * 2 + 0) << mV_extracted.col(k), col;
			outputMeshWireframe.col(i * 2 * posy + j * 2 + 1) << mV_extracted.col(kn), col;
			O_gpu.col(i*posy + j) = mV_extracted.col(k);
			N_gpu.col(i*posy + j) = mNf_extracted.col(i);
		}
	}

	mOutputMeshShader.bind();
	mOutputMeshShader.uploadAttrib("position", O_gpu);
	mOutputMeshShader.uploadAttrib("normal", N_gpu);
	mOutputMeshShader.uploadIndices(F_gpu);
	mOutputMeshFaces = F_gpu.cols();
	mOutputMeshLines = outputMeshWireframe.cols();
	mOutputMeshWireframeShader.bind();
	mOutputMeshWireframeShader.uploadAttrib("position", MatrixXf(outputMeshWireframe.block(0, 0, 3, mOutputMeshLines)));
	mOutputMeshWireframeShader.uploadAttrib("color", MatrixXf(outputMeshWireframe.block(3, 0, 3, mOutputMeshLines)));

	while (!(mLayers[OutputMesh]->checked() && mLayers[OutputMeshWireframe]->checked()))
		keyboardEvent(GLFW_KEY_BACKSLASH, 0, true, 0);
		*/
}

void Parametrizer::ComputeOrientationSingularities()
{
	const MatrixXf &N = hierarchy.mN[0], &Q = hierarchy.mQ[0];
	const MatrixXi &F = hierarchy.mF;
	singularities.clear();
	for (int f = 0; f < F.cols(); ++f) {
		int index = 0;
		int abs_index = 0;
		for (int k = 0; k < 3; ++k) {
			int i = F(k, f), j = F(k == 2 ? 0 : (k + 1), f);
			auto value = compat_orientation_extrinsic_index_4(Q.col(i), N.col(i), Q.col(j), N.col(j));
			index += value.second - value.first;
			abs_index += std::abs(value.second - value.first);
		}
		int index_mod = modulo(index, 4);
		if (index_mod == 1 || index_mod == 3) {
			singularities[f] = index_mod;
		}
	}
	printf("singularities %d...\n", singularities.size());
}