#include "Parametrizer.h"

#include "loader.h"
#include "MergeVertex.h"
#include "subdivide.h"
#include "dedge.h"
#include "AdjacentMatrix.h"
#include "field_math.h"
#include "mcmf.h"
#include <fstream>
#include <queue>
#include <list>
#include <Eigen/Sparse>
#include <GL/glut.h>

#define LOG_OUTPUT

inline double fast_acos(double x) {
	double negate = double(x < 0.0f);
	x = std::abs(x);
	double ret = -0.0187293f;
	ret *= x; ret = ret + 0.0742610f;
	ret *= x; ret = ret - 0.2121144f;
	ret *= x; ret = ret + 1.5707288f;
	ret = ret * std::sqrt(1.0f - x);
	ret = ret - 2.0f * negate * ret;
	return negate * (double)M_PI + ret;
}

void Parametrizer::Load(const char* filename)
{
	load(filename, V, F);
	double maxV[3] = { -1e30, -1e30, -1e30 };
	double minV[3] = { 1e30, 1e30, 1e30 };
	for (int i = 0; i < V.cols(); ++i) {
		for (int j = 0; j < 3; ++j) {
			maxV[j] = std::max(maxV[j], V(j, i));
			minV[j] = std::min(minV[j], V(j, i));
		}
	}
	double scale = std::max(std::max(maxV[0] - minV[0], maxV[1] - minV[1]), maxV[2] - minV[2]) * 0.5;
	for (int i = 0; i < V.cols(); ++i) {
		for (int j = 0; j < 3; ++j) {
			V(j, i) = (V(j, i) - (maxV[j] + minV[j]) * 0.5) / scale;
		}
	}
#ifdef LOG_OUTPUT
	printf("vertices size: %d\n", V.cols());
	printf("faces size: %d\n", F.cols());
#endif

	merge_close(V, F, 1e-6);

}

void Parametrizer::Initialize()
{
	ComputeMeshStatus();

	num_vertices = V.cols() / 4;
	num_faces = num_vertices;
	scale = sqrt(surface_area / num_faces);
	printf("scale %lf\n", scale);
	printf("Compute Direct Graph\n");
	compute_direct_graph(V, F, V2E, E2E, boundary, nonManifold);
	printf("Compute Direct Graph finish\n");

	double target_len = std::min(scale / 2, average_edge_length * 2);
	if (target_len < max_edge_length) {
		subdivide(F, V, V2E, E2E, boundary, nonManifold, target_len);
		compute_direct_graph(V, F, V2E, E2E, boundary, nonManifold);
	}
	printf("Adjacency Matrix\n");
	generate_adjacency_matrix_uniform(F, V2E, E2E, nonManifold, adj);
	printf("Adjacency Matrix finish\n");
	ComputeSmoothNormal();
	ComputeVertexArea();
	triangle_space.resize(F.cols());
	for (int i = 0; i < F.cols(); ++i) {
		Matrix3d p, q;
		p.col(0) = V.col(F(1, i)) - V.col(F(0, i));
		p.col(1) = V.col(F(2, i)) - V.col(F(0, i));
		p.col(2) = Nf.col(i);
		q = p.inverse();
		triangle_space[i].resize(2, 3);
		for (int j = 0; j < 2; ++j) {
			for (int k = 0; k < 3; ++k) {
				triangle_space[i](j, k) = q(j, k);
			}
		}
	}
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
		Vector3d v[3] = { V.col(F(0, f)), V.col(F(1, f)), V.col(F(2, f)) };
		double area = 0.5f * (v[1] - v[0]).cross(v[2] - v[0]).norm();
		surface_area += area;
		for (int i = 0; i < 3; ++i) {
			double len = (v[(i + 1) % 3] - v[i]).norm();
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
	Nf.resize(3, F.cols());
	for (int f = 0; f < F.cols(); ++f) {
		Vector3d v0 = V.col(F(0, f)),
			v1 = V.col(F(1, f)),
			v2 = V.col(F(2, f)),
			n = (v1 - v0).cross(v2 - v0);
		double norm = n.norm();
		if (norm < RCPOVERFLOW) {
			n = Vector3d::UnitX();
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
			N.col(i) = Vector3d::UnitX();
			continue;
		}

		int stop = edge;
		Vector3d normal = Vector3d::Zero();
		do {
			int idx = edge % 3;

			Vector3d d0 = V.col(F((idx + 1) % 3, edge / 3)) - V.col(i);
			Vector3d d1 = V.col(F((idx + 2) % 3, edge / 3)) - V.col(i);
			double angle = fast_acos(d0.dot(d1) / std::sqrt(d0.squaredNorm() * d1.squaredNorm()));

			/* "Computing Vertex Normals from Polygonal Facets"
			by Grit Thuermer and Charles A. Wuethrich, JGT 1998, Vol 3 */
			if (std::isfinite(angle))
				normal += Nf.col(edge / 3) * angle;

			int opp = E2E[edge];
			if (opp == -1)
				break;

			edge = dedge_next_3(opp);
		} while (edge != stop);
		double norm = normal.norm();
		N.col(i) = norm > RCPOVERFLOW ? Vector3d(normal / norm)
			: Vector3d::UnitX();
	}

}


void Parametrizer::ComputeVertexArea() {
	A.resize(V.cols());
	A.setZero();

	for (int i = 0; i < V2E.size(); ++i) {
		int edge = V2E[i], stop = edge;
		if (nonManifold[i] || edge == -1)
			continue;
		double vertex_area = 0;
		do {
			int ep = dedge_prev_3(edge), en = dedge_next_3(edge);

			Vector3d v = V.col(F(edge % 3, edge / 3));
			Vector3d vn = V.col(F(en % 3, en / 3));
			Vector3d vp = V.col(F(ep % 3, ep / 3));

			Vector3d face_center = (v + vp + vn) * (1.0f / 3.0f);
			Vector3d prev = (v + vp) * 0.5f;
			Vector3d next = (v + vn) * 0.5f;

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


void
extract_graph(const Hierarchy &mRes,
std::vector<std::vector<TaggedLink> > &adj_new,
MatrixXd &O_new, MatrixXd &N_new, std::map<uint32_t, uint32_t>& vertex_map,
bool remove_spurious_vertices,
bool remove_unnecessary_edges, int with_scale) {

	double scale = mRes.mScale, inv_scale = 1 / scale;

	auto compat_orientation = compat_orientation_extrinsic_4;
	auto compat_position = compat_position_extrinsic_index_4;

	const MatrixXd &Q = mRes.mQ[0], &O = mRes.mO[0], &N = mRes.mN[0], &V = mRes.mV[0], &S = mRes.mS[0];
	const AdjacentMatrix &adj = mRes.mAdj[0];


	DisjointSets dset(V.cols());
	adj_new.clear();
	adj_new.resize(V.cols());
	typedef std::pair<Edge, double> WeightedEdge;

	std::vector<WeightedEdge> collapse_edge_vec;
	collapse_edge_vec.reserve((uint32_t)(V.cols()*2.5f));
	for (uint32_t i = 0; i < adj.size(); ++i) {
		while (!dset.try_lock(i))
			;

		for (auto& link : adj[i]) {
			uint32_t j = link.id;

			if (j < i)
				continue;

			std::pair<Vector3d, Vector3d> Q_rot = compat_orientation(
				Q.col(i), N.col(i), Q.col(j), N.col(j));
			auto Q_ind = compat_orientation_extrinsic_index_4(
				Q.col(i), N.col(i), Q.col(j), N.col(j));

			double error = 0;
			std::pair<Vector2i, Vector2i> shift;
			if (with_scale == 0)
				shift = compat_position(
				V.col(i), N.col(i), Q_rot.first, O.col(i),
				V.col(j), N.col(j), Q_rot.second, O.col(j),
				scale, scale, inv_scale, inv_scale,
				scale, scale, inv_scale, inv_scale, &error);
			else {
				float scale_x = scale * S(0, i);
				float scale_y = scale * S(1, i);
				float scale_x_1 = scale * S(0, j);
				float scale_y_1 = scale * S(1, j);
				if (Q_ind.first % 2 != Q_ind.second % 2)
					std::swap(scale_x_1, scale_y_1);
				shift = compat_position(
					V.col(i), N.col(i), Q_rot.first, O.col(i),
					V.col(j), N.col(j), Q_rot.second, O.col(j),
					scale_x, scale_y, 1.0f / scale_x, 1.0f / scale_y,
					scale_x_1, scale_y_1, 1.0f / scale_x_1, 1.0f / scale_y_1, &error);
			}

			Vector2i absDiff = (shift.first - shift.second).cwiseAbs();

			if (absDiff.maxCoeff() > 1 || (absDiff == Vector2i(1, 1)))
				continue; /* Ignore longer-distance links and diagonal lines for quads */
			bool collapse = absDiff.sum() == 0;

			if (collapse) {
				collapse_edge_vec.push_back(std::make_pair(std::make_pair(i, j), error));
			}
			else {
				while (!dset.try_lock(j))
					;
				adj_new[i].push_back(j);
				adj_new[j].push_back(i);
				dset.unlock(j);
			}

		}
		dset.unlock(i);

	}

	struct WeightedEdgeComparator {
		bool operator()(const WeightedEdge& e1, const WeightedEdge& e2) const { return e1.second < e2.second; }
	};

	std::stable_sort(collapse_edge_vec.begin(), collapse_edge_vec.end(), WeightedEdgeComparator());

	std::atomic<uint32_t> nConflicts(0), nItem(0);
	std::vector<uint16_t> nCollapses(V.cols(), 0);
	std::set<uint32_t> temp;

	for (int i = 0; i < collapse_edge_vec.size(); ++i) {
		const WeightedEdge &we = collapse_edge_vec[nItem++];
		Edge edge = we.first;

		/* Lock both sets and determine the current representative ID */
		bool ignore_edge = false;
		do {
			if (edge.first > edge.second) {
				std::swap(edge.first, edge.second);
			}
			if (!dset.try_lock(edge.first))
				continue;
			if (!dset.try_lock(edge.second)) {
				dset.unlock(edge.first);
				if (edge.second == edge.first) {
					ignore_edge = true;
					break;
				}
				continue;
			}
			break;
		} while (true);

		if (ignore_edge)
			continue;

		bool contained = false;
		for (auto neighbor : adj_new[edge.first]) {
			if (dset.find(neighbor.id) == edge.second) {
				contained = true;
				break;
			}
		}

		if (contained) {
			dset.unlock(edge.first);
			dset.unlock(edge.second);
			nConflicts++;
			continue;
		}

		temp.clear();
		for (auto neighbor : adj_new[edge.first])
			temp.insert(dset.find(neighbor.id));
		for (auto neighbor : adj_new[edge.second])
			temp.insert(dset.find(neighbor.id));

		uint32_t target_idx = dset.unite_index_locked(edge.first, edge.second);
		adj_new[edge.first].clear();
		adj_new[edge.second].clear();
		adj_new[target_idx].reserve(temp.size());
		for (auto j : temp)
			adj_new[target_idx].push_back(j);
		nCollapses[target_idx] = nCollapses[edge.first] + nCollapses[edge.second] + 1;
		adj_new[edge.first].shrink_to_fit();
		adj_new[edge.second].shrink_to_fit();

		dset.unite_unlock(edge.first, edge.second);
	}

	uint32_t nVertices = 0;
	double avg_collapses = 0;
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
		temp.clear();
		for (auto k : adj_new[i])
			temp.insert(vertex_map[dset.find(k.id)]);
		std::vector<TaggedLink> new_vec;
		new_vec.reserve(temp.size());
		for (auto j : temp)
			new_vec.push_back(TaggedLink(j));
		adj_new[i] = std::move(new_vec);
	}

	if (remove_spurious_vertices) {
		uint32_t removed = 0;
		for (uint32_t i = 0; i<adj_new.size(); ++i) {
			if (nCollapses[i] > avg_collapses / 10)
				continue;

			for (auto neighbor : adj_new[i]) {
				auto &a = adj_new[neighbor.id];
				a.erase(std::remove_if(a.begin(), a.end(), [&](const TaggedLink &v) { return v.id == i; }), a.end());
			}

			adj_new[i].clear();
			++removed;
		}
	}

	O_new.resize(3, nVertices);
	N_new.resize(3, nVertices);
	O_new.setZero();
	N_new.setZero();


	{
		Eigen::VectorXd cluster_weight(nVertices);
		cluster_weight.setZero();
		std::map<uint32_t, uint32_t> vertex_map_new;
		for (int i = 0; i < V.cols(); ++i) {
			auto it = vertex_map.find(dset.find(i));
			if (it == vertex_map.end())
				continue;
			uint32_t j = it->second;
			vertex_map_new[i] = j;

			double weight = std::exp(-(O.col(i) - V.col(i)).squaredNorm() * inv_scale * inv_scale * 9);

			for (uint32_t k = 0; k<3; ++k) {
				O_new.coeffRef(k, j) += O(k, i)*weight;
				N_new.coeffRef(k, j) += N(k, i)*weight;
			}
			cluster_weight[j] += weight;
		}
		std::swap(vertex_map, vertex_map_new);
		
		for (uint32_t i = 0; i<nVertices; ++i) {
			if (cluster_weight[i] == 0) {
				continue;
			}
			O_new.col(i) /= cluster_weight[i];
			N_new.col(i).normalize();
		}

	}
	
	if (remove_unnecessary_edges) {
		bool changed;
		uint32_t nRemoved = 0, nSnapped = 0;
		do {
			changed = false;

			bool changed_inner;
			do {
				changed_inner = false;
				double thresh = 0.3f * scale;

				std::vector<std::tuple<double, uint32_t, uint32_t, uint32_t>> candidates;
				for (uint32_t i_id = 0; i_id<adj_new.size(); ++i_id) {
					auto const &adj_i = adj_new[i_id];
					const Vector3d p_i = O_new.col(i_id);
					for (uint32_t j = 0; j<adj_i.size(); ++j) {
						uint32_t j_id = adj_i[j].id;
						const Vector3d p_j = O_new.col(j_id);
						auto const &adj_j = adj_new[j_id];

						for (uint32_t k = 0; k<adj_j.size(); ++k) {
							uint32_t k_id = adj_j[k].id;
							if (k_id == i_id)
								continue;
							const Vector3d p_k = O_new.col(k_id);
							double a = (p_j - p_k).norm(), b = (p_i - p_j).norm(), c = (p_i - p_k).norm();
							if (a > std::max(b, c)) {
								double s = 0.5f * (a + b + c);
								double height = 2 * std::sqrt(s*(s - a)*(s - b)*(s - c)) / a;
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

					const Vector3d p_i = O_new.col(i), p_j = O_new.col(j), p_k = O_new.col(k);
					double a = (p_j - p_k).norm(), b = (p_i - p_j).norm(), c = (p_i - p_k).norm();
					double s = 0.5f * (a + b + c);
					double height = 2 * std::sqrt(s*(s - a)*(s - b)*(s - c)) / a;
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
						Vector3d n_k = N_new.col(k), n_j = N_new.col(j);
						//Vector3d dp = p_k - p_j, dn = n_k - n_j;
						//double t = dp.dot(p_i-p_j) / dp.dot(dp);
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

				std::vector<std::pair<double, Edge>> candidates;
				for (size_t i = 0; i<adj_new.size(); ++i) {
					auto const &adj_i = adj_new[i];
					const Vector3d p_i = O_new.col(i);
					for (uint32_t j = 0; j<adj_i.size(); ++j) {
						uint32_t j_id = adj_i[j].id;
						const Vector3d p_j = O_new.col(j_id);

						uint32_t nTris = 0;
						double length = 0.0f;
						for (uint32_t k = 0; k<adj_i.size(); ++k) {
							uint32_t k_id = adj_i[k].id;
							if (k_id == j_id)
								continue;
							const Vector3d p_k = O_new.col(k_id);
							if (std::find_if(adj_new[j_id].begin(), adj_new[j_id].end(),
								[k_id](const TaggedLink &l) { return l.id == k_id; }) == adj_new[j_id].end())
								continue;
							nTris++;
							length += (p_k - p_i).norm() + (p_k - p_j).norm();
						}

						if (nTris == 2) {
							double exp_diag = length / 4 * std::sqrt(2.f);
							double diag = (p_i - p_j).norm();
							double score = std::abs((diag - exp_diag) / std::min(diag, exp_diag));
							candidates.push_back(std::make_pair(std::abs(score), std::make_pair(i, j_id)));
						}
					}
				std::sort(candidates.begin(), candidates.end(), [&](
					const std::pair<double, Edge> &v0,
					const std::pair<double, Edge> &v1) { return v0.first < v1.first; });

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

	for (int i = 0; i < O_new.cols(); ++i) {
		Vector3d s, t, p = O_new.col(i);
		coordinate_system(N_new.col(i), s, t);

		std::sort(adj_new[i].begin(), adj_new[i].end(),
			[&](const TaggedLink &j0, const TaggedLink &j1) {
			Vector3d v0 = O_new.col(j0.id) - p, v1 = O_new.col(j1.id) - p;
			return std::atan2(t.dot(v0), s.dot(v0)) > std::atan2(t.dot(v1), s.dot(v1));
		}
		);
	}
}

void extract_faces(std::vector<std::vector<TaggedLink> > &adj, MatrixXd &O,
	MatrixXd &N, MatrixXd &Nf, MatrixXi &F,
	double scale, bool fill_holes) {

	uint32_t nF = 0, nV = O.cols(), nV_old = O.cols();
	F.resize(4, O.cols());

	auto extract_face = [&](uint32_t cur, uint32_t curIdx, size_t targetSize,
		std::vector<std::pair<uint32_t, uint32_t>> &result) {
		uint32_t initial = cur;
		bool success = false;
		result.clear();
		for (;;) {
			if (adj[cur][curIdx].used() ||
				(targetSize > 0 && result.size() + 1 > targetSize))
				break;

			result.push_back(std::make_pair(cur, curIdx));

			uint32_t next = adj[cur][curIdx].id,
				next_rank = adj[next].size(),
				idx = (uint32_t)-1;

			for (uint32_t j = 0; j<next_rank; ++j) {
				if (adj[next][j].id == cur) {
					idx = j; break;
				}
			}

			if (idx == (uint32_t)-1 || next_rank == 1)
				break;

			cur = next;
			curIdx = (idx + 1) % next_rank;
			if (cur == initial) {
				success = targetSize == 0 || result.size() == targetSize;
				break;
			}
		}

		if (success) {
			for (auto kv : result)
				adj[kv.first][kv.second].markUsed();
		}
		else {
			result.clear();
		}
		return success;
	};

	std::vector<std::vector<uint32_t>> irregular_faces;
	auto fill_face = [&](std::vector<std::pair<uint32_t, uint32_t>> &verts) -> std::vector<uint32_t> {
		std::vector<uint32_t> irregular;
		while (verts.size() > 2) {
			if (verts.size() == (size_t)4) {
				uint32_t idx = nF++;
				if (nF > F.cols())
					F.conservativeResize(F.rows(), F.cols() * 2);
				for (int i = 0; i < 4; ++i)
					F(i, idx) = verts[i].first;
				break;
			}
			else if (verts.size() >(size_t) 4 + 1 || 4 == 3) {
				double best_score = std::numeric_limits<double>::infinity();
				uint32_t best_idx = (uint32_t)-1;

				for (uint32_t i = 0; i<verts.size(); ++i) {
					double score = 0.f;
					for (int k = 0; k < 4; ++k) {
						Vector3d v0 = O.col(verts[(i + k) % verts.size()].first);
						Vector3d v1 = O.col(verts[(i + k + 1) % verts.size()].first);
						Vector3d v2 = O.col(verts[(i + k + 2) % verts.size()].first);
						Vector3d d0 = (v0 - v1).normalized();
						Vector3d d1 = (v2 - v1).normalized();
						double angle = std::acos(d0.dot(d1)) * 180.0f / M_PI;
						score += std::abs(angle - (4 == 4 ? 90 : 60));
					}

					if (score < best_score) {
						best_score = score;
						best_idx = i;
					}
				}
				uint32_t idx = nF++;
				if (nF > F.cols())
					F.conservativeResize(F.rows(), F.cols() * 2);

				for (int i = 0; i < 4; ++i) {
					uint32_t &j = verts[(best_idx + i) % verts.size()].first;
					F(i, idx) = j;
					if (i != 0 && (int)i != 4 - 1)
						j = (uint32_t)-1;
				}
				verts.erase(std::remove_if(verts.begin(), verts.end(),
					[](const std::pair<uint32_t, uint32_t> &v) { return v.first == (uint32_t)-1; }), verts.end());
			}
			else {
				Vector3d centroid = Vector3d::Zero();
				Vector3d centroid_normal = Vector3d::Zero();
				for (uint32_t k = 0; k<verts.size(); ++k) {
					centroid += O.col(verts[k].first);
					centroid_normal += N.col(verts[k].first);
				}
				uint32_t idx_centroid = nV++;
				if (nV > O.cols()) {
					O.conservativeResize(O.rows(), O.cols() * 2);
					N.conservativeResize(O.rows(), O.cols());
				}
				O.col(idx_centroid) = centroid / verts.size();
				N.col(idx_centroid) = centroid_normal.normalized();

				for (uint32_t i = 0; i<verts.size(); ++i) {
					uint32_t idx = nF++;
					if (nF > F.cols())
						F.conservativeResize(F.rows(), F.cols() * 2);

					F.col(idx) = Vector4i(
						verts[i].first,
						verts[(i + 1) % verts.size()].first,
						idx_centroid,
						idx_centroid
						);
					irregular.push_back(idx);
				}
				break;
			}
		}
		return irregular;
	};

	VectorXi stats(10);
	stats.setZero();
	uint32_t nFaces = 0, nHoles = 0;
	std::vector<std::pair<uint32_t, uint32_t>> result;
	for (uint32_t _deg = 3; _deg <= 8; _deg++) {
		uint32_t deg = _deg;
		if ((deg == 3 || deg == 4))
			deg = 7 - deg;

		for (uint32_t i = 0; i<nV_old; ++i) {
			for (uint32_t j = 0; j<adj[i].size(); ++j) {
				if (!extract_face(i, j, _deg, result))
					continue;
				stats[result.size()]++;
				std::vector<uint32_t> irregular = fill_face(result);
				if (!irregular.empty())
					irregular_faces.push_back(std::move(irregular));
				nFaces++;
			}
		}
	}

	if (fill_holes) {
		for (uint32_t i = 0; i<nV_old; ++i) {
			for (uint32_t j = 0; j<adj[i].size(); ++j) {
				if (!adj[i][j].used()) {
					uint32_t j_id = adj[i][j].id;
					bool found = false;
					for (uint32_t k = 0; k<adj[j_id].size(); ++k) {
						if (adj[j_id][k].id == i) {
							found = true;
							if (adj[j_id][k].used()) {
								adj[i][j].flag |= 2;
								adj[j_id][k].flag |= 2;
							}
							break;
						}
					}
				}
			}
		}

		uint32_t linksLeft = 0;
		for (uint32_t i = 0; i<nV_old; ++i) {
			adj[i].erase(std::remove_if(adj[i].begin(), adj[i].end(),
				[](const TaggedLink &l) { return (l.flag & 2) == 0; }), adj[i].end());
			linksLeft += adj[i].size();
		}

		for (uint32_t i = 0; i<nV_old; ++i) {
			for (uint32_t j = 0; j<adj[i].size(); ++j) {
				if (!extract_face(i, j, 0, result))
					continue;
				if (result.size() >= 7) {
					continue;
				}
				if (result.size() >= (size_t)stats.size()) {
					uint32_t oldSize = stats.size();
					stats.conservativeResize(result.size() + 1);
					stats.tail(stats.size() - oldSize).setZero();
				}
				stats[result.size()]++;
				std::vector<uint32_t> irregular = fill_face(result);
				if (!irregular.empty())
					irregular_faces.push_back(std::move(irregular));
				nHoles++;
			}
		}
	}

	F.conservativeResize(4, nF);
	N.conservativeResize(3, nV);
	O.conservativeResize(3, nV);

	Nf.resize(3, F.cols());
	Nf.setZero();
	for (uint32_t i = 0; i < F.cols(); ++i) {
		Vector3d centroid = Vector3d::Zero(), avgNormal = Vector3d::Zero();
		for (int j = 0; j<F.rows(); ++j) {
			uint32_t k = F(j, i);
			centroid += O.col(k);
			avgNormal += N.col(k);
		}
		centroid /= F.rows();
		Matrix3d cov = Matrix3d::Zero();
		for (int j = 0; j<F.rows(); ++j) {
			uint32_t k = F(j, i);
			cov += (O.col(k) - centroid) * (O.col(k) - centroid).transpose();
		}
		Vector3d n = cov.jacobiSvd(Eigen::ComputeFullU).matrixU().col(2).normalized();
		Nf.col(i) = n * signum(avgNormal.dot(n));
	}

	for (auto f : irregular_faces) {
		Vector3d centroid = Vector3d::Zero(), avgNormal = Vector3d::Zero();
		for (uint32_t i = 0; i<f.size(); ++i) {
			uint32_t k = F(0, f[i]);
			centroid += O.col(k);
			avgNormal += N.col(k);
		}
		centroid /= f.size();
		Matrix3d cov = Matrix3d::Zero();
		for (uint32_t i = 0; i<f.size(); ++i) {
			uint32_t k = F(0, f[i]);
			cov += (O.col(k) - centroid) * (O.col(k) - centroid).transpose();
		}
		Vector3d n = cov.jacobiSvd(Eigen::ComputeFullU).matrixU().col(2).normalized();
		n *= signum(avgNormal.dot(n));
		for (uint32_t i = 0; i<f.size(); ++i)
			Nf.col(f[i]) = n;
	}
}

void write_obj(const std::string &filename, const MatrixXi &F,
	const MatrixXd &V, const MatrixXd &N, const MatrixXd &Nf,
	const MatrixXd &UV, const MatrixXd &C) {
	std::ofstream os(filename);
	if (os.fail())
		throw std::runtime_error("Unable to open OBJ file \"" + filename + "\"!");
	if (N.size() > 0 && Nf.size() > 0)
		throw std::runtime_error("Please specify either face or vertex normals but not both!");

	for (uint32_t i = 0; i<V.cols(); ++i)
		os << "v " << V(0, i) << " " << V(1, i) << " " << V(2, i) << std::endl;

	for (uint32_t i = 0; i<N.cols(); ++i)
		os << "vn " << N(0, i) << " " << N(1, i) << " " << N(2, i) << std::endl;

	for (uint32_t i = 0; i<Nf.cols(); ++i)
		os << "vn " << Nf(0, i) << " " << Nf(1, i) << " " << Nf(2, i) << std::endl;

	for (uint32_t i = 0; i<UV.cols(); ++i)
		os << "vt " << UV(0, i) << " " << UV(1, i) << std::endl;

	/* Check for irregular faces */
	std::map<uint32_t, std::pair<uint32_t, std::map<uint32_t, uint32_t>>> irregular;
	size_t nIrregular = 0;

	for (uint32_t f = 0; f<F.cols(); ++f) {
		if (F.rows() == 4) {
			if (F(2, f) == F(3, f)) {
				nIrregular++;
				auto &value = irregular[F(2, f)];
				value.first = f;
				value.second[F(0, f)] = F(1, f);
				continue;
			}
		}
		os << "f ";
		for (uint32_t j = 0; j<F.rows(); ++j) {
			uint32_t idx = F(j, f);
			idx += 1;
			os << idx;
			if (Nf.size() > 0)
				idx = f + 1;
			os << "//" << idx << " ";
		}
		os << std::endl;
	}

	printf("irregular faces: %d\n", irregular.size());
	for (auto item : irregular) {
		auto face = item.second;
		uint32_t v = face.second.begin()->first, first = v, i = 0;
		os << "f ";
		while (true) {
			uint32_t idx = v + 1;
			os << idx;
			if (Nf.size() > 0)
				idx = face.first + 1;
			os << "//" << idx << " ";

			v = face.second[v];
			if (v == first || ++i == face.second.size())
				break;
		}
		os << std::endl;
	}

}

void Parametrizer::ComputeOrientationSingularities()
{
	MatrixXd &N = hierarchy.mN[0], &Q = hierarchy.mQ[0];
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
			if (index >= 4 || index < 0) {
				Q.col(F(0, f)) = -Q.col(F(0, f));
			}
			singularities[f] = index_mod;
		}
	}
	printf("singularities %d...\n", singularities.size());
}

void Parametrizer::ComputePositionSingularities(int with_scale)
{
	const MatrixXd &V = hierarchy.mV[0], &N = hierarchy.mN[0], &Q = hierarchy.mQ[0], &O = hierarchy.mO[0];
	const MatrixXi &F = hierarchy.mF;

	pos_sing.clear();
	pos_rank.resize(F.rows(), F.cols());
	pos_index.resize(6, F.cols());
	for (int f = 0; f < F.cols(); ++f) {
		if (singularities.count(f))
			continue;

		Vector2i index = Vector2i::Zero();
		uint32_t i0 = F(0, f), i1 = F(1, f), i2 = F(2, f);

		Vector3d q[3] = { Q.col(i0).normalized(), Q.col(i1).normalized(), Q.col(i2).normalized() };
		Vector3d n[3] = { N.col(i0), N.col(i1), N.col(i2) };
		Vector3d o[3] = { O.col(i0), O.col(i1), O.col(i2) };
		Vector3d v[3] = { V.col(i0), V.col(i1), V.col(i2) };

		int best[3];
		double best_dp = -std::numeric_limits<double>::infinity();
		for (int i = 0; i<4; ++i) {
			Vector3d v0 = rotate90_by(q[0], n[0], i);
			for (int j = 0; j<4; ++j) {
				Vector3d v1 = rotate90_by(q[1], n[1], j);
				for (int k = 0; k<4; ++k) {
					Vector3d v2 = rotate90_by(q[2], n[2], k);
					double dp = std::min(std::min(v0.dot(v1), v1.dot(v2)), v2.dot(v0));
					if (dp > best_dp) {
						best_dp = dp;
						best[0] = i; best[1] = j; best[2] = k;
					}
				}
			}
		}
		pos_rank(0, f) = best[0];
		pos_rank(1, f) = best[1];
		pos_rank(2, f) = best[2];
		for (int k = 0; k<3; ++k)
			q[k] = rotate90_by(q[k], n[k], best[k]);

		for (int k = 0; k<3; ++k) {
			int kn = k == 2 ? 0 : (k + 1);
			double scale_x = hierarchy.mScale, scale_y = hierarchy.mScale, scale_x_1 = hierarchy.mScale, scale_y_1 = hierarchy.mScale;
			if (with_scale) {
				scale_x *= hierarchy.mS[0](0, F(k, f));
				scale_y *= hierarchy.mS[0](1, F(k, f));
				scale_x_1 *= hierarchy.mS[0](0, F(kn, f));
				scale_y_1 *= hierarchy.mS[0](1, F(kn, f));
				if (best[k] % 2 != 0)
					std::swap(scale_x, scale_y);
//				if (best[k] >= 2)
//					scale_x = -scale_x, scale_y = -scale_y;
				if (best[kn] % 2 != 0)
					std::swap(scale_x_1, scale_y_1);
//				if (best[kn] >= 2)
//					scale_x_1 = -scale_x_1, scale_y_1 = -scale_y_1;
			}
			double inv_scale_x = 1.0 / scale_x, inv_scale_y = 1.0 / scale_y, inv_scale_x_1 = 1.0 / scale_x_1, inv_scale_y_1 = 1.0 / scale_y_1;
			std::pair<Vector2i, Vector2i> value =
				compat_position_extrinsic_index_4(
				v[k], n[k], q[k], o[k],
				v[kn], n[kn], q[kn], o[kn],
				scale_x, scale_y, inv_scale_x, inv_scale_y,
				scale_x_1, scale_y_1, inv_scale_x_1, inv_scale_y_1, nullptr);
			auto diff = value.first - value.second;
			index += diff;
			pos_index(k * 2, f) = diff[0];
			pos_index(k * 2 + 1, f) = diff[1];
		}

		if (index != Vector2i::Zero()) {
			pos_sing[f] = rshift90(index, best[0]);
		}
	}

	printf("position singularities %d...\n", pos_sing.size());
}

void Parametrizer::ExtractMesh(int with_scale) {
	printf("extract_graph\n");
	std::map<uint32_t, uint32_t> vertex_map;
	extract_graph(hierarchy, adj_extracted,
		mV_extracted, mN_extracted, vertex_map, true, true, with_scale);
	MatrixXi F_extr;
	MatrixXd Nf_extr;

	auto adj_new = adj_extracted;
	extract_faces(adj_new, mV_extracted, mN_extracted, Nf_extr, F_extr,
		hierarchy.mScale, true);

	write_obj("result.obj", F_extr, mV_extracted, MatrixXd(), Nf_extr, MatrixXd(), MatrixXd());
	printf("finish writing...\n");
	std::map<uint32_t, std::pair<uint32_t, std::map<uint32_t, uint32_t>>> irregular;
	size_t nIrregular = 0;

	qF.clear();
	for (uint32_t f = 0; f<F_extr.cols(); ++f) {
		if (F_extr.rows() == 4) {
			if (F_extr(2, f) == F_extr(3, f)) {
				nIrregular++;
				auto &value = irregular[F_extr(2, f)];
				value.first = f;
				value.second[F_extr(0, f)] = F_extr(1, f);
				continue;
			}
		}
		Vector4i p = F_extr.col(f);
		VectorXi px(4);
		for (int i = 0; i < 4; ++i) {
			px(i) = p(i);
		}
		qF.push_back(px);
	}

	for (auto item : irregular) {
		auto face = item.second;
		uint32_t v = face.second.begin()->first, first = v, i = 0;
		std::vector<uint32_t> p;
		while (true) {
			p.push_back(v);

			v = face.second[v];
			if (v == first || ++i == face.second.size())
				break;
		}
		VectorXi px(p.size());
		for (int i = 0; i < p.size(); ++i) {
			px(i) = p[i];
		}
		qF.push_back(px);
	}

	qE.clear();
	qVF.resize(mV_extracted.size());
	for (int i = 0; i < qF.size(); ++i) {
		std::vector<int> cw_edge_ids;
		std::vector<int> ccw_edge_ids;
		for (int j = 0; j < qF[i].rows(); ++j) {
			int v0 = qF[i](j);
			int v1 = qF[i]((j + 1) % qF[i].rows());
			qVF[v0].push_back(std::make_pair(i, j));
			Edge e(v0, v1);
			int idx = qE.size();
			if (edge_idmap.count(e) == 0) {
				edge_idmap[e] = qE.size();
				qE.push_back(e);
				triangle_edge_pair.push_back(std::set<int>());
			}
			else {
				idx = edge_idmap[e];
			}
			cw_edge_ids.push_back(idx);
			idx = qE.size();
			e = Edge(v1, v0);
			if (edge_idmap.count(e) == 0) {
				edge_idmap[e] = qE.size();
				qE.push_back(e);
				triangle_edge_pair.push_back(std::set<int>());
			}
			else {
				idx = edge_idmap[e];
			}
			ccw_edge_ids.push_back(idx);
		}

		if (qF[i].rows() == 3) {
			int v0 = qF[i](0);
			int v1 = qF[i](1);
			int v2 = qF[i](2);
			double angle1 = (mV_extracted.col(v1) - mV_extracted.col(v0)).normalized().dot(
				(mV_extracted.col(v2) - mV_extracted.col(v0)).normalized());
			double angle2 = (mV_extracted.col(v2) - mV_extracted.col(v1)).normalized().dot(
				(mV_extracted.col(v0) - mV_extracted.col(v1)).normalized());
			double angle3 = (mV_extracted.col(v0) - mV_extracted.col(v2)).normalized().dot(
				(mV_extracted.col(v1) - mV_extracted.col(v2)).normalized());

			int edge1, edge2;
			if (angle1 > angle2 && angle1 > angle3) {
				edge1 = edge_idmap[Edge(v0, v1)];
				edge2 = edge_idmap[Edge(v0, v2)];
				triangle_edge_pair[edge1].insert(edge2);
				triangle_edge_pair[edge2].insert(edge1);
			}
			else if (angle2 > angle1 && angle2 > angle3) {
				edge1 = edge_idmap[Edge(v1, v0)];
				edge2 = edge_idmap[Edge(v1, v2)];
				triangle_edge_pair[edge1].insert(edge2);
				triangle_edge_pair[edge2].insert(edge1);
			}
			else {
				edge1 = edge_idmap[Edge(v2, v0)];
				edge2 = edge_idmap[Edge(v2, v1)];
				triangle_edge_pair[edge1].insert(edge2);
				triangle_edge_pair[edge2].insert(edge1);
			}
		}
		
		if (qF[i].rows() == 5) {
			double max_dot = -2.0f;
			int index = -1, opposite_index = -1;
			for (int j = 0; j < 5; ++j) {
				int va = qF[i]((j + 4) % qF[i].rows());
				int v0 = qF[i](j);
				int v1 = qF[i]((j + 1) % qF[i].rows());
				Vector3d p = mV_extracted.col(v0) - mV_extracted.col(va);
				Vector3d q = mV_extracted.col(v1) - mV_extracted.col(v0);
				double dot = p.normalized().dot(q.normalized());
				if (dot > max_dot) {
					max_dot = dot;
					index = j;
				}
			}
			int candidate1 = qF[i]((index + 2) % 5);
			int candidate2 = qF[i]((index + 3) % 5);
			int current_v = qF[i](index);
			int next_v = qF[i]((index + 1) % 5);
			int prev_v = qF[i]((index + 4) % 5);
			Vector3d p = (mV_extracted.col(candidate1) - mV_extracted.col(current_v)).normalized();
			Vector3d q = (mV_extracted.col(candidate2) - mV_extracted.col(current_v)).normalized();
			Vector3d t = (mV_extracted.col(next_v) - mV_extracted.col(current_v)).normalized();
			Edge e(current_v, candidate2);
			if (edge_idmap.count(e) == 0) {
				edge_idmap[e] = qE.size();
				qE.push_back(e);
				triangle_edge_pair.push_back(std::set<int>());
			}
			e = Edge(e.second, e.first);
			if (edge_idmap.count(e) == 0) {
				edge_idmap[e] = qE.size();
				qE.push_back(e);
				triangle_edge_pair.push_back(std::set<int>());
			}
			int edge1 = edge_idmap[Edge(e.second, next_v)];
			int edge2 = edge_idmap[e];
			triangle_edge_pair[edge1].insert(edge2);
			triangle_edge_pair[edge2].insert(edge1);
		}
	}

	qVE.resize(mV_extracted.cols());
	qEV.resize(mV_extracted.cols());
	for (int e = 0; e < qE.size(); ++e) {
		qVE[qE[e].first].push_back(e);
		qEV[qE[e].second].push_back(e);
	}

	qRE.resize(qE.size(), -1);
	for (int i = 0; i < qE.size(); ++i) {
		Edge e(qE[i].second, qE[i].first);
		auto it = edge_idmap.find(e);
		if (it != edge_idmap.end()) {
			qRE[i] = it->second;
		}
	}

	qEE.resize(qE.size(), -1);
	for (int i = 0; i < qE.size(); ++i) {
		Vector3d p = mV_extracted.col(qE[i].second) - mV_extracted.col(qE[i].first);
		p = p.cross(Vector3d(mN_extracted.col(qE[i].second))).normalized();

		int reverse_edge_id = qRE[i];
		double max_dot = 0.0f;
		for (auto e : qVE[qE[i].second]) {
			bool found = false;
			if (e == reverse_edge_id)
				found = true;
			if (!found) {
				Vector3d q = Vector3d(mV_extracted.col(qE[e].second) - mV_extracted.col(qE[e].first)).cross(
					Vector3d(mN_extracted.col(qE[i].second))).normalized();
				double dot = p.dot(q.normalized());
				if (dot > max_dot) {
					max_dot = dot;
					qEE[i] = e;
				}
			}
		}
	}

	std::map<int, int> extracted_singularities;
	for (int i = 0; i < mV_extracted.cols(); ++i) {
		if (qVE[i].size() % 4 != 0) {
			extracted_singularities[i] = 1;
		}
	}

	for (int i = 0; i < qF.size(); ++i) {
		if (qF[i].rows() != 4) {
			for (int j = 0; j < qF[i].rows(); ++j) {
				extracted_singularities[qF[i](j)] = 3;
			}
		}
	}
	
	for (auto& f : singularities) {
		int face = f.first;
		bool find_singularity = false;
		for (int j = 0; j < 3; ++j) {
			int v = hierarchy.mF(j, face);
			int qv = vertex_map[v];
			if (extracted_singularities.count(qv)) {
				vertex_singularities[qv] = f.second;
				find_singularity = true;
				break;
			}
		}
		if (find_singularity == false) {
			for (int k = 0; k < 3; ++k)
				vertex_singularities[vertex_map[hierarchy.mF(k, face)]] = f.second;
		}
	}
}

void Parametrizer::LoopFace(int mode)
{
	if (mode == 0 || mode == 2) {
		sin_graph.resize(mV_extracted.cols());
		q.reserve(1000000);
		double angle_thres = 45.0;
		for (auto& s : vertex_singularities) {
			// loop the edges
			std::vector<Vector3d> directions;
			std::vector<int> edge_ids;
			int q_index = q.size();
			for (int i = 0; i < qVE[s.first].size(); ++i) {
				int edge_id = qVE[s.first][i];
				Vector3d dir = mV_extracted.col(qE[edge_id].second) - mV_extracted.col(qE[edge_id].first);
				dir.normalize();
				bool found = false;
				for (int j = 0; j < directions.size(); ++j) {
					if (//fabs((directions[j].dot(dir))) > cos(angle_thres / 180.0*3.141592654) ||
						triangle_edge_pair[edge_id].count(edge_ids[j])) {
						found = true;
						break;
					}
				}
				if (!found) {
					ExpandInfo info;
					info.current_v = s.first;
					info.singularity = s.first;
					info.edge_id = edge_id;
					info.step = 0;
					info.prev = -1;
					q.push_back(info);
					directions.push_back(dir);
					edge_ids.push_back(edge_id);
				}
			}
			sin_graph[s.first][s.first] = std::make_pair(-1, q_index);
		}
		front_index = 0;
		singularity_entry.resize(mV_extracted.cols());
	}
	if (mode == 0)
		return;

	auto GenEdgeStrip = [&](int qid, std::vector<int>& edges, int reverse) {
		while (true) {
			if (reverse == 1)
				edges.push_back(qRE[q[qid].edge_id]);
			else
				edges.push_back(q[qid].edge_id);
			qid = q[qid].prev;
			if (qid == -1)
				break;
		}
	};
	auto TryInsert = [&](int v1, int e1, int v2, int e2, std::vector<int>& edges) {
		bool found = singularity_entry[v1].count(std::make_pair(v2, e1)) || singularity_entry[v2].count(std::make_pair(v1, e2));
		if (!found) {
			for (auto& e : triangle_edge_pair[e1]) {
				if (singularity_entry[v1].count(std::make_pair(v2, e))) {
					found = true;
					break;
				}
			}
			for (auto& e : triangle_edge_pair[e2]) {
				if (singularity_entry[v2].count(std::make_pair(v1, e))) {
					found = true;
					break;
				}
			}
		}

		if (!found) {
			singularity_entry[v1].insert(std::make_pair(v2, e1));
			for (auto& e : triangle_edge_pair[e1]) {
				singularity_entry[v1].insert(std::make_pair(v2, e));
			}
			singularity_entry[v2].insert(std::make_pair(v1, e2));
			for (auto& e : triangle_edge_pair[e2]) {
				singularity_entry[v2].insert(std::make_pair(v1, e));
			}
			edge_strips.push_back(std::move(edges));
		}
	};
	while (front_index < q.size()) {
		auto& info = q[front_index];
		int next_v = qE[info.edge_id].second;
		int next_edge_id = qEE[info.edge_id];
		bool find_terminal = false;
		if (!sin_graph[next_v].empty()) {
			// generate a strip
			for (auto& p : sin_graph[next_v]) {
				auto& next_info = q[p.second.second];
				int next_edge = p.second.first;
				if (next_edge == -1 || qEE[next_edge] == qRE[info.edge_id] || qEE[info.edge_id] == qRE[next_edge]) {
					std::vector<int> edge_id1, edge_id2;
					GenEdgeStrip(front_index, edge_id1, 0);
					GenEdgeStrip(p.second.second, edge_id2, 1);
					std::reverse(edge_id1.begin(), edge_id1.end());
					edge_id1.insert(edge_id1.end(), edge_id2.begin() + 1, edge_id2.end());

					int v1 = qE[edge_id1.front()].first, e1 = edge_id1.front();
					int v2 = qE[edge_id1.back()].second, e2 = qRE[edge_id1.back()];
					TryInsert(v1, e1, v2, e2, edge_id1);
					find_terminal = true;
				}
				else
				if (next_info.step <= 1) {
					// principal direction is edge_id
					Vector3d dir = (mV_extracted.col(qE[info.edge_id].second) - mV_extracted.col(qE[info.edge_id].first)).normalized();
					double max_dot = 2.0f;
					int edge_entry = 0;
					for (int i = 0; i < qVE[p.first].size(); ++i) {
						int e = qVE[p.first][i];
						Vector3d next_dir = (mV_extracted.col(qE[e].second) - mV_extracted.col(qE[e].first)).normalized();
						double dot = dir.dot(next_dir);
						if (dot < max_dot) {
							max_dot = dot;
							edge_entry = e;
						}
					}
					std::vector<int> edge_id1, edge_id2;
					GenEdgeStrip(front_index, edge_id1, 0);
					GenEdgeStrip(p.second.second, edge_id2, 1);
					std::reverse(edge_id1.begin(), edge_id1.end());

					edge_id1.insert(edge_id1.end(), edge_id2.begin() + 1, edge_id2.end());

					int v1 = qE[edge_id1.front()].first, e1 = edge_id1.front();
					int v2 = qE[edge_id1.back()].second, e2 = edge_entry;
					TryInsert(v1, e1, v2, e2, edge_id1);
					find_terminal = true;
				}
			}
		}
		if (sin_graph[next_v].count(info.singularity) != 0) {
			find_terminal = true;
		}
		else {
			sin_graph[next_v][info.singularity] = std::make_pair(info.edge_id, q.size());
		}
		if (!find_terminal && next_edge_id != -1) {
			ExpandInfo new_info;
			new_info.current_v = next_v;
			new_info.singularity = info.singularity;
			new_info.edge_id = next_edge_id;
			new_info.step = info.step + 1;
			new_info.prev = front_index;
			q.push_back(new_info);
		}
		front_index += 1;
		if (mode == 1)
			break;
	}
}

void Parametrizer::EstimateScale() {
	auto& mF = hierarchy.mF;
	auto& mQ = hierarchy.mQ[0];
	auto& mN = hierarchy.mN[0];
	auto& mV = hierarchy.mV[0];
	FS.resize(2, mF.cols());
	FQ.resize(3, mF.cols());
	for (int i = 0; i < mF.cols(); ++i) {
		const Vector3d& n = Nf.col(i);
		const Vector3d &q_1 = mQ.col(mF(0, i)), &q_2 = mQ.col(mF(1, i)), &q_3 = mQ.col(mF(2, i));
		const Vector3d &n_1 = mN.col(mF(0, i)), &n_2 = mN.col(mF(1, i)), &n_3 = mN.col(mF(2, i));
		Vector3d q_1n = rotate_vector_into_plane(q_1, n_1, n);
		Vector3d q_2n = rotate_vector_into_plane(q_2, n_2, n);
		Vector3d q_3n = rotate_vector_into_plane(q_3, n_3, n);

		auto p = compat_orientation_extrinsic_4(q_1n, n, q_2n, n);
		Vector3d q = (p.first + p.second).normalized();
		p = compat_orientation_extrinsic_4(q, n, q_3n, n);
		q = (p.first * 2 + p.second);
		q = q - n * q.dot(n);
		FQ.col(i) = q.normalized();
	}
	for (int i = 0; i < mF.cols(); ++i) {
		double step = hierarchy.mScale * 1.f;

		const Vector3d &n = Nf.col(i);
		Vector3d p = (mV.col(mF(0, i)) + mV.col(mF(1, i)) + mV.col(mF(2, i))) * (1.0 / 3.0);
		Vector3d q_x = FQ.col(i), q_y = n.cross(q_x);
		Vector3d q_xl = -q_x, q_xr = q_x;
		Vector3d q_yl = -q_y, q_yr = q_y;
		Vector3d q_yl_unfold = q_y, q_yr_unfold = q_y, q_xl_unfold = q_x, q_xr_unfold = q_x;
		int f;
		double tx, ty, len;

		f = i; len = step;
		TravelField(p, q_xl, len, f, hierarchy.mE2E, mV, mF, Nf, FQ, mQ, mN, triangle_space, &tx, &ty, &q_yl_unfold);

		f = i; len = step;
		TravelField(p, q_xr, len, f, hierarchy.mE2E, mV, mF, Nf, FQ, mQ, mN, triangle_space, &tx, &ty, &q_yr_unfold);

		f = i; len = step;
		TravelField(p, q_yl, len, f, hierarchy.mE2E, mV, mF, Nf, FQ, mQ, mN, triangle_space, &tx, &ty, &q_xl_unfold);

		f = i; len = step;
		TravelField(p, q_yr, len, f, hierarchy.mE2E, mV, mF, Nf, FQ, mQ, mN, triangle_space, &tx, &ty, &q_xr_unfold);
		double dSx = (q_yr_unfold - q_yl_unfold).dot(q_x) / (2.0f * step);
		double dSy = (q_xr_unfold - q_xl_unfold).dot(q_y) / (2.0f * step);
		FS.col(i) = Vector2d(dSx, dSy);
	}
	
	std::vector<double> areas(mV.cols(), 0.0);
	for (int i = 0; i < mF.cols(); ++i) {
		Vector3d p1 = mV.col(mF(1, i)) - mV.col(mF(0, i));
		Vector3d p2 = mV.col(mF(2, i)) - mV.col(mF(0, i));
		double area = p1.cross(p2).norm();
		for (int j = 0; j < 3; ++j) {
			auto index = compat_orientation_extrinsic_index_4(FQ.col(i), Nf.col(i), mQ.col(mF(j, i)), mN.col(mF(j, i)));
			double scaleX = FS.col(i).x(), scaleY = FS.col(i).y();
			if (index.first != index.second % 2) {
				std::swap(scaleX, scaleY);
			}
			if (index.second >= 2) {
				scaleX = -scaleX;
				scaleY = -scaleY;
			}
			hierarchy.mK[0].col(mF(j, i)) += area * Vector2d(scaleX, scaleY);
			areas[mF(j, i)] += area;
		}
	}
	for (int i = 0; i < mV.cols(); ++i) {
		if (areas[i] != 0)
			hierarchy.mK[0].col(i) /= areas[i];
	}
	for (int l = 0; l< hierarchy.mK.size() - 1; ++l)  {
		const MatrixXd &K = hierarchy.mK[l];
		MatrixXd &K_next = hierarchy.mK[l + 1];
		auto& toUpper = hierarchy.mToUpper[l];
		for (int i = 0; i < toUpper.cols(); ++i) {
			Vector2i upper = toUpper.col(i);
			Vector2d k0 = K.col(upper[0]);

			if (upper[1] != -1) {
				Vector2d k1 = K.col(upper[1]);
				k0 = 0.5 * (k0 + k1);
			}

			K_next.col(i) = k0;
		}
	}
}

void Parametrizer::SaveToFile(FILE* fp) {
	Save(fp, vertex_singularities);
	Save(fp, singularities);
	Save(fp, V);
	Save(fp, N);
	Save(fp, Nf);
	Save(fp, FQ);
	Save(fp, FS);
	Save(fp, F);
	Save(fp, V2E);
	Save(fp, E2E);
	Save(fp, boundary);
	Save(fp, nonManifold);
	Save(fp, adj);
	Save(fp, triangle_space);
	hierarchy.SaveToFile(fp);
	Save(fp, surface_area);
	Save(fp, scale);
	Save(fp, average_edge_length);
	Save(fp, max_edge_length);
	Save(fp, A);
	Save(fp, num_vertices);
	Save(fp, adj_extracted);
	Save(fp, mF_extracted);
	Save(fp, mV_extracted);
	Save(fp, mN_extracted);
	Save(fp, mNf_extracted);
	Save(fp, qF);
	Save(fp, qVF);
	Save(fp, edge_idmap);
	Save(fp, qE);
	Save(fp, qVE);
	Save(fp, qEV);
	Save(fp, qEE);
	Save(fp, qRE);
	Save(fp, sin_graph);
	Save(fp, triangle_edge_pair);
}

void Parametrizer::LoadFromFile(FILE* fp) {
	Read(fp, vertex_singularities);
	Read(fp, singularities);
	Read(fp, V);
	Read(fp, N);
	Read(fp, Nf);
	Read(fp, FQ);
	Read(fp, FS);
	Read(fp, F);
	Read(fp, V2E);
	Read(fp, E2E);
	Read(fp, boundary);
	Read(fp, nonManifold);
	Read(fp, adj);
	Read(fp, triangle_space);
	hierarchy.LoadFromFile(fp);
	Read(fp, surface_area);
	Read(fp, scale);
	Read(fp, average_edge_length);
	Read(fp, max_edge_length);
	Read(fp, A);
	Read(fp, num_vertices);
	Read(fp, adj_extracted);
	Read(fp, mF_extracted);
	Read(fp, mV_extracted);
	Read(fp, mN_extracted);
	Read(fp, mNf_extracted);
	Read(fp, qF);
	Read(fp, qVF);
	Read(fp, edge_idmap);
	Read(fp, qE);
	Read(fp, qVE);
	Read(fp, qEV);
	Read(fp, qEE);
	Read(fp, qRE);
	Read(fp, sin_graph);
	Read(fp, triangle_edge_pair);
}

int get_parents(std::vector<std::pair<int, int> >& parents, int j) {
	if (j == parents[j].first)
		return j;
	int k = get_parents(parents, parents[j].first);
	parents[j].second = (parents[j].second + parents[parents[j].first].second) % 4;
	parents[j].first = k;
	return k;
}

int get_parents_orient(std::vector<std::pair<int, int> >& parents, int j) {
	if (j == parents[j].first)
		return parents[j].second;
	return (parents[j].second + get_parents_orient(parents, parents[j].first)) % 4;
}

void Parametrizer::MergeVertices(int v) {
	std::set<DEdge> different_pairs;
	{
		std::vector<std::map<int, std::pair<int, Vector2i> > > links(shrink_parents.size());
		for (auto& s : singular_patches)
			if (singular_patches.count(get_parents(shrink_parents, s)) == 0) {
				printf("fail %d %d!\n", s, get_parents(shrink_parents, s));
				system("pause");
			}
		for (auto& e : Es) {
			int v1 = get_parents(shrink_parents, e.first.first);
			int v2 = get_parents(shrink_parents, e.first.second);
			if (singular_patches.count(v1) || singular_patches.count(v2))
				continue;
			int orient = get_parents_orient(shrink_parents, e.first.first);
			auto diff = rshift90(e.second.second, orient);
			// vi - vo = vi - e.first
			int diff_orient = (-orient + e.second.first + get_parents_orient(shrink_parents, e.first.second) + 4) % 4;
			if (links[v1].count(v2)) {
				if (links[v1][v2].second != diff || links[v1][v2].first != diff_orient)
					links[v1][v2].second = Vector2i(10000, 10000);
			}
			else {
				links[v1][v2] = std::make_pair(diff_orient, diff);
			}
		}
		Vector2i entries[] = { Vector2i(1, 0), Vector2i(0, 1), Vector2i(-1, 0), Vector2i(0, -1) };
		auto shrink_parents_buf = shrink_parents;
		for (int i = 0; i < links.size(); ++i) {
			std::list<int> l[4];
			for (int j = 0; j < 4; ++j) {
				for (auto& e : links[i]) {
					if (e.second.second == entries[j]) {
						l[j].push_back(get_parents(shrink_parents, e.first));
					}
				}
			}
			for (int j = 0; j < 4; ++j) {
				for (auto& p1 : l[j]) {
					for (int k = j + 1; k < 4; ++k) {
						for (auto& p2 : l[k]) {
							different_pairs.insert(DEdge(p1, p2));
						}
					}
					different_pairs.insert(DEdge(p1, i));
				}
			}
		}
	}
	std::map<int, std::pair<double, Vector2i> > links;
	for (auto& e : Es) {
		int v1 = shrink_compact_indices[get_parents(shrink_parents, e.first.first)];
		int v2 = get_parents(shrink_parents, e.first.second);
		if (v1 != v)
			continue;
		if (singular_patches.count(v1) || singular_patches.count(v2))
			continue;
		int orient = get_parents_orient(shrink_parents, e.first.first);
		auto diff = rshift90(e.second.second, orient);
		// vi - vo = vi - e.first
		int diff_orient = (-orient + e.second.first + get_parents_orient(shrink_parents, e.first.second) + 4) % 4;
		if (links.count(v2)) {
			if (links[v2].second != diff || links[v2].first != diff_orient)
				links[v2].second = Vector2i(10000, 10000);
		}
		else {
			links[v2] = std::make_pair(diff_orient, diff);
		}
	}
	Vector2i entries[] = { Vector2i(1, 0), Vector2i(0, 1), Vector2i(-1, 0), Vector2i(0, -1) };
	auto shrink_parents_buf = shrink_parents;
	for (int j = 0; j < 4; ++j) {
		std::list<std::pair<int, int> > l;
		for (auto& e : links) {
			if (e.second.second == entries[j]) {
				l.push_back(std::make_pair(e.first, e.second.first));
			}
		}
		printf("l.size %d\n", l.size());
		if (l.empty())
			continue;
		auto& v = l.front();
		for (auto& e : l) {
			int v0 = get_parents(shrink_parents, v.first);
			int v1 = get_parents(shrink_parents, e.first);
			if (v0 != v1) {
				if (different_pairs.count(DEdge(v0, v1))) {
					printf("%d %d %d %d\n", v.first, e.first, v0, v1);
					printf("continue...\n");
					continue;
				}
				if (singular_patches.count(v1) == 0 && (shrink_ranks[v0] > shrink_ranks[v1] || singular_patches.count(v0) != 0)) {
					shrink_ranks[v0] += shrink_ranks[v1];
					int diff = (get_parents_orient(shrink_parents, v0) + v.second - e.second - get_parents_orient(shrink_parents, v1) + 8) % 4;
					shrink_parents[v1] = std::make_pair(v0, diff);
				}
				else {
					shrink_ranks[v1] += shrink_ranks[v0];
					int diff = (get_parents_orient(shrink_parents, v1) + e.second - v.second - get_parents_orient(shrink_parents, v0) + 8) % 4;
					shrink_parents[v0] = std::make_pair(v1, diff);
				}
			}
			else {
				if (v.second != e.second) {
					singular_patches.insert(v0);
					printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
				}
			}
		}
	}
	compact_num_v = 0;
	for (int i = 0; i < shrink_parents.size(); ++i) {
		if (shrink_parents[i].first == i)
			shrink_compact_indices[i] = compact_num_v++;
	}
}

void Parametrizer::UpdateMesh()
{
	mE_extracted.clear();
	mE_extracted2.clear();
	mO_extracted.clear();
	mO_extracted.resize(compact_num_v, Vector3d(0, 0, 0));
	std::vector<double> shrink_counter(compact_num_v, 0);
	for (int i = 0; i < O_compact.size(); ++i) {
		int p = shrink_compact_indices[get_parents(shrink_parents, i)];
		shrink_counter[p] += counter[i];
		mO_extracted[p] += O_compact[i] * counter[i];
	}
	for (int i = 0; i < mO_extracted.size(); ++i) {
		mO_extracted[i] /= shrink_counter[i];
	}

	singular_patches_buf.clear();
	singular_e_buf.clear();
	for (auto& s : singular_patches) {
		singular_patches_buf.insert(shrink_compact_indices[get_parents(shrink_parents, s)]);
	}

	std::set<DEdge> edges, edges2;
	
	for (auto& e : Es) {
		int v1 = shrink_compact_indices[get_parents(shrink_parents, e.first.first)];
		int v2 = shrink_compact_indices[get_parents(shrink_parents, e.first.second)];
		if (v1 != v2) {
			if (abs(e.second.second[0]) + abs(e.second.second[1]) > 1)
				edges2.insert(DEdge(v1, v2));
			else
				edges.insert(DEdge(v1, v2));
		}
	}
	/*
	for (int i = 0; i < edge_values.size(); ++i) {
		int v1 = shrink_compact_indices[get_parents(shrink_parents, colors[edge_values[i].x])];
		int v2 = shrink_compact_indices[get_parents(shrink_parents, colors[edge_values[i].y])];
		if (v1 != v2) {
			if (abs(edge_diff[i][0]) + abs(edge_diff[i][1]) > 1)
				edges2.insert(DEdge(v1, v2));
			else
				edges.insert(DEdge(v1, v2));
		}
	}
	*/
	for (auto& e : singular_e) {
		auto diff = edge_diff[edge_ids[e]];
		int v1 = shrink_compact_indices[get_parents(shrink_parents, e.x)];
		int v2 = shrink_compact_indices[get_parents(shrink_parents, e.y)];
		singular_e_buf.insert(DEdge(v1, v2));
	}
	mE_extracted.insert(mE_extracted.end(), edges.begin(), edges.end());
	mE_extracted2.insert(mE_extracted2.end(), edges2.begin(), edges2.end());
	auto& F = hierarchy.mF;
	auto& Q = hierarchy.mQ[0];
	auto& N = hierarchy.mN[0];
	
	for (int i = 0; i < F.cols(); ++i) {
		if (singularities.count(i))
			continue;
		int shift_edges = 0;
		{
			for (int j = 0; j < 3; ++j) {
				int v1 = F(j, i);
				int v2 = F((j + 1) % 3, i);
				int v3 = F((j + 2) % 3, i);
				auto diff1 = edge_diff[edge_ids[DEdge(v1, v2)]];
				auto diff2 = edge_diff[edge_ids[DEdge(v1, v3)]];
				auto rank1 = compat_orientation_extrinsic_index_4(Q.col(v1), N.col(v1), Q.col(v2), N.col(v2));
				auto rank2 = compat_orientation_extrinsic_index_4(Q.col(v1), N.col(v1), Q.col(v3), N.col(v3));
				if (v1 > v2)
					diff1 = rshift90(-diff1, (rank1.first + 4 - rank1.second) % 4);
				if (v1 > v3)
					diff2 = rshift90(-diff2, (rank2.first + 4 - rank2.second) % 4);
				if (diff1[0] * diff2[1] - diff1[1] * diff2[0] < 0) {
					shift_edges += 1;
//					printf("<%d %d> => <%d %d>\n", diff1[0], diff1[1], diff2[0], diff2[1]);
				}
			}
			if (shift_edges > 0) {
//				printf("################# shift edges %d #####################\n", i);
			}
		}
		int v1 = shrink_compact_indices[get_parents(shrink_parents, colors[F(0, i)])];
		int v2 = shrink_compact_indices[get_parents(shrink_parents, colors[F(1, i)])];
		int v3 = shrink_compact_indices[get_parents(shrink_parents, colors[F(2, i)])];
		if (v1 != v2 && v2 != v3 && v3 != v1 && shift_edges == 1) {
//			printf("shifts %d\n", shift_edges);
//			mF_extracted2.push_back(Vector3i(v1, v2, v3));
		}
	}
}

void Parametrizer::ComputeIndexMap(int with_scale)
{
	// compute indices
	auto& V = hierarchy.mV[0];
	auto& F = hierarchy.mF;
	auto& adj = hierarchy.mAdj[0];
	auto& Q = hierarchy.mQ[0];
	auto& N = hierarchy.mN[0];
	auto& O = hierarchy.mO[0];
	auto &S = hierarchy.mS[0];

	edge_diff.clear();
	edge_ids.clear();
	edge_neighbors.reserve(F.cols() * 2);
	for (int i = 0; i < F.cols(); ++i) {
		if (singularities.count(i)) {
			for (int j = 0; j < 3; ++j) {
				int v1 = F(j, i);
				int v2 = F((j + 1) % 3, i);
				singular_edges.insert(DEdge(v1, v2));
			}
			continue;
		}
		for (int j = 0; j < 3; ++j) {
			int k0 = (j + 2) % 3, k1 = j, k2 = (j + 1) % 3;
			int v0 = F(k0, i);
			int v1 = F(k1, i);
			int v2 = F(k2, i);
			DEdge e1(v0, v1), e2(v1, v2);
			Vector2i diff1, diff2;
			int rank1, rank2;
			if (v0 > v1) {
				rank1 = pos_rank(k1, i);
				diff1 = rshift90(Vector2i(-pos_index(k0 * 2, i), -pos_index(k0 * 2 + 1, i)), rank1);
			}
			else {
				rank1 = pos_rank(k0, i);
				diff1 = rshift90(Vector2i(pos_index(k0 * 2, i), pos_index(k0 * 2 + 1, i)), rank1);
			}
			if (v1 > v2) {
				rank2 = pos_rank(k2, i);
				diff2 = rshift90(Vector2i(-pos_index(k1 * 2, i), -pos_index(k1 * 2 + 1, i)), rank2);
			}
			else {
				rank2 = pos_rank(k1, i);
				diff2 = rshift90(Vector2i(pos_index(k1 * 2, i), pos_index(k1 * 2 + 1, i)), rank2);
			}

			int eID1, eID2;
			auto it = edge_ids.find(e1);
			if (it == edge_ids.end()) {
				edge_ids[e1] = edge_values.size();
				eID1 = edge_values.size();
				edge_values.push_back(e1);
				edge_neighbors.push_back(std::list<int>());
				edge_neighbors.push_back(std::list<int>());
				edge_diff.push_back(diff1);
			}
			else {
				eID1 = it->second;
			}
			eID1 = eID1 * 2 + (v0 > v1);
			it = edge_ids.find(e2);
			if (it == edge_ids.end()) {
				edge_ids[e2] = edge_values.size();
				eID2 = edge_values.size();
				edge_values.push_back(e2);
				edge_neighbors.push_back(std::list<int>());
				edge_neighbors.push_back(std::list<int>());
				edge_diff.push_back(diff2);
			}
			else {
				eID2 = it->second;
			}
			eID2 = eID2 * 2 + (v1 > v2);
			edge_neighbors[eID1].push_back(eID2);
			edge_neighbors[eID2].push_front(eID1);
		}
	}

	/*
	int count1 = 0, count2 = 0;
	for (int i = 0; i < F.cols(); ++i) {
		if (singularities.count(i))
			continue;
		int v1 = F(0, i);
		int v2 = F(1, i);
		int v3 = F(2, i);
		std::vector<int> ts = { v1, v2, v3 };
		std::sort(ts.begin(), ts.end());
		if (ts[0] == 189 && ts[1] == 190 && ts[2] == 9920) {
			ts = ts;
		}
		auto diff1 = edge_diff[edge_ids[DEdge(v1, v2)]];
		auto diff2 = edge_diff[edge_ids[DEdge(v1, v3)]];
		auto index1 = compat_orientation_extrinsic_index_4(Q.col(v1), N.col(v1), Q.col(v2), N.col(v2));
		auto index2 = compat_orientation_extrinsic_index_4(Q.col(v1), N.col(v1), Q.col(v3), N.col(v3));
		int rank1 = (index1.first - index1.second + 4) % 4;
		int rank2 = (index2.first - index2.second + 4) % 4;
		if (v1 > v2)
			diff1 = rshift90(-diff1, rank1);
		if (v1 > v3)
			diff2 = rshift90(-diff2, rank2);
		int t = diff1[0] * diff2[1] - diff1[1] * diff2[0];
		if (t > 0)
			count1 += 1;
		if (t < 0) {
			count2 += 1;
			if (pos_sing.count(i))
				printf("1\n");
			else
				printf("0\n");
		}
	}
	printf("counter %d %d\n", count1, count2);
	system("pause");
	*/
	
	auto RankDiff = [&](int e1, int e2) {
		int v0 = edge_values[e1].x;
		int v1 = edge_values[e2].x;
		auto t = compat_orientation_extrinsic_index_4(Q.col(v0), N.col(v0),
			Q.col(v1), N.col(v1));
		return (t.second + 4 - t.first) % 4;
	};

	auto TransformRank = [&](Vector2i& d1, int e1, int e2) {
		auto diff = RankDiff(e1, e2);
		return rshift90(d1, diff);
	};

	auto TransformRankV = [&](const Vector2i& d1, int v1, int v2) {
		auto t = compat_orientation_extrinsic_index_4(Q.col(v1), N.col(v1),
			Q.col(v2), N.col(v2));
		return rshift90(d1, (t.second + 4 - t.first) % 4);
	};

	printf("compute hash map\n");
	int t1 = GetTickCount(), t2 = t1;
	std::vector<int> edge_hash_map(edge_values.size() * 2, -2);
	for (int f = 0; f < F.cols(); ++f) {
		if (singularities.count(f)) {
			for (int i = 0; i < 3; ++i) {
				DEdge e(F(i, f), F((i + 1) % 3, f));
				int edgeID = edge_ids[e];
				edge_hash_map[edgeID * 2] = f;
				edge_hash_map[edgeID * 2 + 1] = f;
			}
		}
	}
	auto CheckManifold = [&](const Vector2i& potential_diff, int eID) {
		if (edge_neighbors[eID].empty())
			return true;
		int id1 = eID / 2;
		int id2 = edge_neighbors[eID].front() / 2;
		auto& e1 = edge_values[id1];
		auto& e2 = edge_values[id2];
		int intersect;
		if (e1.x == e2.x || e1.x == e2.y)
			intersect = e1.x;
		if (e1.y == e2.x || e1.y == e2.y)
			intersect = e1.y;
		auto diff1 = potential_diff;//edge_diff[id1];
		auto diff2 = edge_diff[id2];
		if (intersect == e1.y) {
			auto index = compat_orientation_extrinsic_index_4(Q.col(e1.x), N.col(e1.x), Q.col(e1.y), N.col(e1.y));
			diff1 = rshift90(-diff1, (index.second - index.first + 4) % 4);
		}
		if (intersect == e2.y) {
			auto index = compat_orientation_extrinsic_index_4(Q.col(e2.x), N.col(e2.x), Q.col(e2.y), N.col(e2.y));
			diff2 = rshift90(-diff2, (index.second - index.first + 4) % 4);
		}
		return (diff1[0] * diff2[1] - diff1[1] * diff2[0] >= 0);
	};
	auto CheckSanity = [&](int sing){
		printf("Check Sanity\n");
		int pos_sings = 0;
		for (int ff = 0; ff < F.cols(); ++ff) {
			//			if (singularities.count(ff))
			//				continue;
			for (int j = 0; j < 3; ++j) {
				int v0 = F(j, ff);
				int v1 = F((j + 1) % 3, ff);
				int v2 = F((j + 2) % 3, ff);
				auto diff1 = edge_diff[edge_ids[DEdge(v0, v1)]];
				auto diff2 = edge_diff[edge_ids[DEdge(v1, v2)]];
				auto diff3 = edge_diff[edge_ids[DEdge(v2, v0)]];
				auto index1 = compat_orientation_extrinsic_index_4(Q.col(v0), N.col(v0), Q.col(v1), N.col(v1));
				auto index2 = compat_orientation_extrinsic_index_4(Q.col(v0), N.col(v0), Q.col(v2), N.col(v2));
				int rank1 = (index1.first - index1.second + 4) % 4;
				int rank2 = (index2.first - index2.second + 4) % 4;
				if (v1 < v0) {
					diff1 = rshift90(-diff1, rank1);
				}
				if (v2 < v1) {
					diff2 = rshift90(-diff2, rank2);
				}
				else {
					diff2 = rshift90(diff2, rank1);
				}
				if (v0 < v2) {
					diff3 = -diff3;
				}
				else {
					diff3 = rshift90(diff3, rank2);
				}
				auto diff = diff1 + diff2 + diff3;

				if (singularities.count(ff)) {
//					printf("%d (%d): <%d %d  %d %d  %d %d>\n", ff, singularities[ff],
						//					shrink_compact_indices[get_parents(shrink_parents, colors[F(0, ff)])],
						//					shrink_compact_indices[get_parents(shrink_parents, colors[F(1, ff)])],
						//					shrink_compact_indices[get_parents(shrink_parents, colors[F(2, ff)])],
//						diff1[0], diff1[1], diff2[0], diff2[1], diff3[0], diff3[1]);
				}
				if (diff != Vector2i::Zero() && singularities.count(ff) == 0) {
					if (pos_sing.count(ff) == 0) {
						printf("additional %d !\n", ff);
						printf("%d %d %d\n", F(0, ff), F(1, ff), F(2, ff));
					}
					pos_sings += 1;
				}
			}
		}
		int shift_edges = 0;
		for (int f = 0; f < F.cols(); ++f) {
			if (pos_sing.count(f))
				continue;
			for (int j = 0; j < 3; ++j) {
				int v1 = F(j, f);
				int v2 = F((j + 1) % 3, f);
				int v3 = F((j + 2) % 3, f);
				auto diff1 = edge_diff[edge_ids[DEdge(v1, v2)]];
				auto diff2 = edge_diff[edge_ids[DEdge(v1, v3)]];
				auto rank1 = compat_orientation_extrinsic_index_4(Q.col(v1), N.col(v1), Q.col(v2), N.col(v2));
				auto rank2 = compat_orientation_extrinsic_index_4(Q.col(v1), N.col(v1), Q.col(v3), N.col(v3));
				if (v1 > v2)
					diff1 = rshift90(-diff1, (rank1.first + 4 - rank1.second) % 4);
				if (v1 > v3)
					diff2 = rshift90(-diff2, (rank2.first + 4 - rank2.second) % 4);
				if (diff1[0] * diff2[1] - diff1[1] * diff2[0] < 0)
					shift_edges += 1;
			}
		}
		printf("\nsingularity: %d    shift edges: %d    pos sing: %d\n", sing, shift_edges, pos_sings);
	};
	auto CollapseVertex = [&](int v1, int v2) {
		Vector2i diff1 = edge_diff[edge_ids[DEdge(v1, v2)]];
		auto rank1 = compat_orientation_extrinsic_index_4(Q.col(v1), N.col(v1), Q.col(v2), N.col(v2));
		if (v1 > v2)
			diff1 = rshift90(-diff1, (rank1.first + 4 - rank1.second) % 4);
		std::map<int, int> parents;
		std::queue<int> q;
		q.push(v1);
		parents[v1] = 0;
		while (!q.empty()) {
			int v1 = q.front();
			q.pop();
			for (auto& link : adj[v1]) {
				auto& diff_v1 = edge_diff[edge_ids[DEdge(v1, link.id)]];
				if (diff_v1[0] == 0 && diff_v1[1] == 0) {
					if (!parents.count(link.id)) {
						auto index = compat_orientation_extrinsic_index_4(Q.col(v1), N.col(v1), Q.col(link.id), N.col(link.id));
						q.push(link.id);
						parents[link.id] = (parents[v1] + index.second - index.first + 4) % 4;
					}
				}
			}
		}
		for (auto v1 : parents) {
			for (auto& link : adj[v1.first]) {
				Vector2i& diff_v1 = edge_diff[edge_ids[DEdge(v1.first, link.id)]];
				auto rank2 = compat_orientation_extrinsic_index_4(Q.col(v1.first), N.col(v1.first), Q.col(link.id), N.col(link.id));
				auto diff2 = rshift90(diff1, v1.second);
				if (link.id < v1.first)
					diff_v1 = rshift90(-diff_v1, (rank2.first + 4 - rank2.second) % 4);
				diff_v1 -= diff2;
				if (link.id < v1.first)
					diff_v1 = rshift90(-diff_v1, (rank2.second - rank2.first + 4) % 4);
			}
		}

	};
	CheckSanity(0);
	system("pause");
	printf("Build Integer Constraints...\n");
	BuildIntegerConstraints();
	CheckSanity(0);
	system("pause");
	/*
	int sing = 0;
	for (int f = 0; f < F.cols(); ++f) {
		if (pos_sing.count(f) == 0)
			continue;
		int v0 = F(0, f);
		int v1 = F(1, f);
		int v2 = F(2, f);
		DEdge e0(v0, v1), e1(v1, v2), e2(v2, v0);
		int sign[] = { 1, 1, 1 };
		if (v1 < v0)
			sign[0] = -1;
		if (v2 < v1)
			sign[1] = -1;
		if (v0 < v2)
			sign[2] = -1;
		int eID0 = edge_ids[e0], eID1 = edge_ids[e1], eID2 = edge_ids[e2];
		Vector2i diff = (sign[0] * edge_diff[eID0] + sign[1] * TransformRank(edge_diff[eID1], eID1, eID0) + sign[2] * TransformRank(edge_diff[eID2], eID2, eID0));
		if (diff == Vector2i::Zero())// && singularities.count(f) == 0 || (abs(diff[0]) <= 1 && abs(diff[1]) <= 1 && singularities.count(f)))
			continue;
		int min_len = 0x7fffffff;
		int min_sign = 1;
		int edge_id = -1;
		int prev_id = -1;

		struct EdgeInfo
		{
			int current_edge;
			int prev_edge;
			int direction;
			Vector2i current_diff;
		};
		std::vector<EdgeInfo> q;
		std::vector<int> edge_hash = edge_hash_map;
		auto TraceBack = [&](int index) {
			while (index != -1) {
				edge_diff[q[index].current_edge] += q[index].current_diff;
				index = q[index].prev_edge;
			}
		};
		bool found = false;
		for (int i = 0; i < 3; ++i) {
			int v0 = F(i, f);
			int v1 = F((i + 1) % 3, f);
			int eID = edge_ids[DEdge(v0, v1)];
			if (i == 0)
				prev_id = eID;
			auto new_diff = -sign[i] * TransformRank(diff, prev_id, eID);
			auto potential_diff = new_diff + edge_diff[eID];
			if (abs(potential_diff[0]) < 2 || abs(potential_diff[1]) < 2) {
				if (CheckManifold(potential_diff, eID * 2 + (v0 > v1))) {
					EdgeInfo info;
					info.current_edge = eID;
					info.prev_edge = -1;
					info.current_diff = new_diff;
					info.direction = eID * 2 + (v0 > v1);
					q.push_back(info);
					if (edge_hash[info.direction] >= 0) {
						int f = edge_hash[info.direction];
						int v3 = F(0, f) + F(1, f) + F(2, f) - v0 - v1;
						Vector2i diff;
						if (v3 > edge_values[eID].x) {
							diff = edge_diff[eID];
						}
						else {
							auto index = compat_orientation_extrinsic_index_4(Q.col(v3), N.col(v3), Q.col(edge_values[eID].x), N.col(edge_values[eID].x));
							diff = rshift90(-edge_diff[eID], (index.second - index.first + 4) % 4);
						}
						int sign = potential_diff[0] * diff[1] - potential_diff[1] * diff[0];
						if (v1 < v0 && sign >= 0 || v1 > v0 && sign <= 0) {
							found = true;
							TraceBack(q.size() - 1);
							break;
						}
						else {
							q.pop_back();
						}
					}
					edge_hash[eID * 2] = -1;
					edge_hash[eID * 2 + 1] = -1;
				}
			}
		}
		int front = 0;
		while (front < q.size()) {
			if (found)
				break;
			auto info = q[front++];
			int current_eID = info.direction;
			if (current_eID % 2 == 0)
				current_eID += 1;
			else
				current_eID -= 1;
			edge_diff[info.current_edge] += info.current_diff;
			for (auto& n : edge_neighbors[current_eID]) {
				if (edge_hash[n] == -1)
					continue;
				int next_edge = n / 2;
				auto& e1 = edge_values[info.current_edge];
				auto& e2 = edge_values[next_edge];
				int rank_diff = RankDiff(info.current_edge, next_edge);
				Vector2i new_diff;
				if (e1.x == e2.x || e1.y == e2.y) {
					new_diff = rshift90(info.current_diff, rank_diff);
				}
				else {
					new_diff = -rshift90(info.current_diff, rank_diff);
				}
				auto potential_diff = edge_diff[next_edge] + new_diff;
				if (abs(potential_diff[0]) < 2 && abs(potential_diff[1]) < 2) {
					if (CheckManifold(potential_diff, n)) {
						EdgeInfo new_info;
						new_info.current_edge = next_edge;
						new_info.prev_edge = front - 1;
						new_info.current_diff = new_diff;
						new_info.direction = n;
						q.push_back(new_info);
						if (edge_hash[n] >= 0) {
							edge_diff[next_edge] += new_diff;
							int f = edge_hash[n];
							int v0 = F(0, f);
							int v1 = F(1, f);
							int v2 = F(2, f);
							auto index1 = compat_orientation_extrinsic_index_4(Q.col(v0), N.col(v0), Q.col(v1), N.col(v1));
							auto index2 = compat_orientation_extrinsic_index_4(Q.col(v0), N.col(v0), Q.col(v2), N.col(v2));
							auto diff1 = edge_diff[edge_ids[DEdge(v0, v1)]];
							auto diff2 = edge_diff[edge_ids[DEdge(v1, v2)]];
							auto diff3 = edge_diff[edge_ids[DEdge(v2, v0)]];
							int rank1 = (index1.first - index1.second + 4) % 4;
							int rank2 = (index2.first - index2.second + 4) % 4;
							if (v1 < v0) {
								diff1 = rshift90(-diff1, rank1);
							}
							if (v2 < v1) {
								diff2 = rshift90(-diff2, rank2);
							}
							else {
								diff2 = rshift90(diff2, rank1);
							}
							if (v0 < v2) {
								diff3 = -diff3;
							}
							else {
								diff3 = rshift90(diff3, rank2);
							}


							bool flag = false;
							int k1 = abs(diff1[0]) + abs(diff1[1]);
							int k2 = abs(diff2[0]) + abs(diff2[1]);
							int k3 = abs(diff3[0]) + abs(diff3[1]);
							if (k1 == 0 || k2 == 0 || k3 == 0)
								flag = true;

							int v3 = F(0, f) + F(1, f) + F(2, f) - v0 - v1;
							Vector2i diff;
							if (v3 > edge_values[next_edge].x) {
								diff = edge_diff[next_edge];
							}
							else {
								auto index = compat_orientation_extrinsic_index_4(Q.col(v3), N.col(v3), Q.col(edge_values[next_edge].x), N.col(edge_values[next_edge].x));
								diff = rshift90(-edge_diff[next_edge], (index.second - index.first + 4) % 4);
							}
							int sign = potential_diff[0] * diff[1] - potential_diff[1] * diff[0];
							if (!(v1 < v0 && sign >= 0 || v1 > v0 && sign <= 0)) {
								flag = false;
							}
							edge_diff[next_edge] -= new_diff;
							if (flag) {
								found = true;
								TraceBack(q.size() - 1);
								break;
							}
							else {
								q.pop_back();
							}
						}
						edge_hash[next_edge * 2] = -1;
						edge_hash[next_edge * 2 + 1] = -1;
					}
				}
			}
			edge_diff[info.current_edge] -= info.current_diff;
		}
		if (!found) {
			printf("No solution!\n");
//			system("pause");
		}
		sing += 1;
	}
	printf("sing...\n");
	system("pause");

	for (int i = 0; i < F.cols(); ++i) {
		if (singularities.count(i))
			continue;
		int shift_edges = 0;
		for (int j = 0; j < 3; ++j) {
			int v1 = F(j, i);
			int v2 = F((j + 1) % 3, i);
			int v3 = F((j + 2) % 3, i);
			auto diff1 = edge_diff[edge_ids[DEdge(v1, v2)]];
			auto diff2 = edge_diff[edge_ids[DEdge(v1, v3)]];
			auto rank1 = compat_orientation_extrinsic_index_4(Q.col(v1), N.col(v1), Q.col(v2), N.col(v2));
			auto rank2 = compat_orientation_extrinsic_index_4(Q.col(v1), N.col(v1), Q.col(v3), N.col(v3));
			if (v1 > v2)
				diff1 = rshift90(-diff1, (rank1.first + 4 - rank1.second) % 4);
			if (v1 > v3)
				diff2 = rshift90(-diff2, (rank2.first + 4 - rank2.second) % 4);
			if (diff1[0] * diff2[1] - diff1[1] * diff2[0] < 0) {
				shift_edges += 1;
			}
		}
		if (shift_edges > 0) {
			int v1 = F(0, i);
			int v2 = F(1, i);
			int v3 = F(2, i);
			auto diff1 = edge_diff[edge_ids[DEdge(v1, v2)]];
			auto diff2 = edge_diff[edge_ids[DEdge(v2, v3)]];
			auto diff3 = edge_diff[edge_ids[DEdge(v3, v1)]];
			int t1 = abs(diff1[0]) + abs(diff1[1]);
			int t2 = abs(diff2[0]) + abs(diff2[1]);
			int t3 = abs(diff3[0]) + abs(diff3[1]);
			if (t2 <= t1 && t2 <= t3) {
				v1 = v2, v2 = v3;
				diff1 = diff2;
			}
			else if (t3 <= t1 && t3 <= t2) {
				v2 = v1, v1 = v3;
				diff1 = diff3;
			}
			CollapseVertex(v1, v2);
		}
	}
//	CheckSanity(0);
	int t2 = GetTickCount();
	printf("Trace Map use %lf seconds\n", (t2 - t1) * 1e-3);
	*/
	std::vector<Eigen::Triplet<double> > lhsTriplets;
	t1 = GetTickCount();
	lhsTriplets.reserve(F.cols() * 6);
	std::vector<std::map<int, double> > entries(V.cols() * 2);
	std::vector<double> r_entries(V.cols() * 2);
	for (auto& info : edge_ids) {
		int v1 = info.first.x;
		int v2 = info.first.y;
		Vector3d q_1 = Q.col(v1);
		Vector3d q_2 = Q.col(v2);
		Vector3d n_1 = N.col(v1);
		Vector3d n_2 = N.col(v2);
		Vector3d q_1_y = n_1.cross(q_1);
		Vector3d q_2_y = n_2.cross(q_2);
		Vector3d weights[] = { q_2, q_2_y, -q_1, -q_1_y };
		auto index = compat_orientation_extrinsic_index_4(q_1, n_1, q_2, n_2);
		double s_x1 = S(0, v1), s_y1 = S(1, v1);
		double s_x2 = S(0, v2), s_y2 = S(1, v2);
		int rank_diff = (index.second + 4 - index.first) % 4;
		Vector3d qd_x = 0.5 * (rotate90_by(q_2, n_2, rank_diff) + q_1);
		Vector3d qd_y = 0.5 * (rotate90_by(q_2_y, n_2, rank_diff) + q_1_y);
		double scale_x = (with_scale ? 0.5 * (s_x1 + s_x2) : 1) * hierarchy.mScale;
		double scale_y = (with_scale ? 0.5 * (s_y1 + s_y2) : 1) * hierarchy.mScale;
		Vector2i diff = edge_diff[info.second];
		Vector3d C = diff[0] * scale_x * qd_x + diff[1] * scale_y * qd_y + V.col(v1) - V.col(v2);
		int vid[] = { v2 * 2, v2 * 2 + 1, v1 * 2, v1 * 2 + 1 };
		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j) {
				auto it = entries[vid[i]].find(vid[j]);
				if (it == entries[vid[i]].end()) {
					entries[vid[i]][vid[j]] = weights[i].dot(weights[j]);
				}
				else {
					entries[vid[i]][vid[j]] += weights[i].dot(weights[j]);
				}
			}
			r_entries[vid[i]] += weights[i].dot(C);
		}
	}

	Eigen::SparseMatrix<double> A(V.cols() * 2, V.cols() * 2);
	VectorXd rhs(V.cols() * 2);
	rhs.setZero();
	for (int i = 0; i < entries.size(); ++i) {
		rhs(i) = r_entries[i];
		for (auto& rec : entries[i]) {
			lhsTriplets.push_back(Eigen::Triplet<double>(i, rec.first, rec.second));
		}
	}
	A.setFromTriplets(lhsTriplets.begin(), lhsTriplets.end());
	Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
	solver.analyzePattern(A);
	solver.factorize(A);

	VectorXd result = solver.solve(rhs);
	for (int i = 0; i < O.cols(); ++i) {
		Vector3d q = Q.col(i);
		Vector3d n = N.col(i);
		Vector3d q_y = n.cross(q);
		O.col(i) = V.col(i) + q * result[i * 2] + q_y * result[i * 2 + 1];
	}
	
	t2 = GetTickCount();
	printf("Use %lf seconds.\n", (t2 - t1) * 1e-3);

	std::vector<std::pair<int, int> > parents(V.cols());
	std::vector<int> ranks(V.cols(), 1);
	for (int i = 0; i < parents.size(); ++i) {
		parents[i] = std::make_pair(i, 0);
	}

	int edge_num = 0;
	for (int i = 0; i < edge_diff.size(); ++i) {
		if (edge_diff[i] == Vector2i::Zero()) {
			int vv0 = edge_values[i].x;
			int vv1 = edge_values[i].y;
			int v0 = get_parents(parents, edge_values[i].x);
			int v1 = get_parents(parents, edge_values[i].y);
			auto index = compat_orientation_extrinsic_index_4(Q.col(edge_values[i].x),
				N.col(edge_values[i].x), Q.col(edge_values[i].y), N.col(edge_values[i].y));
			int diff = (index.second - index.first + 4) % 4;

			if (v0 == v1) {
				int t = get_parents_orient(parents, vv0) - get_parents_orient(parents, vv1);
				if (diff != (4 + get_parents_orient(parents, vv0) - get_parents_orient(parents, vv1)) % 4) {
				}
				continue;
			}
			if (ranks[v0] > ranks[v1]) {
				ranks[v0] += ranks[v1];
				diff = (8 - diff - get_parents_orient(parents, edge_values[i].y) + get_parents_orient(parents, edge_values[i].x)) % 4;
				parents[v1] = std::make_pair(v0, diff);
			}
			else {
				ranks[v1] += ranks[v0];
				diff = (4 - get_parents_orient(parents, edge_values[i].x) + diff + get_parents_orient(parents, edge_values[i].y)) % 4;
				parents[v0] = std::make_pair(v1, diff);
			}
		}
	}

	printf("quad edge num: %d\n", edge_values.size());
	std::vector<int> compact_parents(parents.size(), -1);
	int num_v;
	int count = 0;
	while (true) {
		std::vector<int> singularity_pos(V.cols(), 0);
		for (auto& f : singularities) {
			for (int i = 0; i < 3; ++i) {
				singularity_pos[get_parents(parents, F(i, f.first))] = 1;
			}
		}
		for (auto& p : cuts) {
			singularity_pos[get_parents(parents, p.x)] = 1;
			singularity_pos[get_parents(parents, p.y)] = 1;
		}

		count += 1;
		bool update = false;
		num_v = 0;
		for (int i = 0; i < parents.size(); ++i) {
			if (parents[i].first == i) {
				compact_parents[i] = num_v++;
			}
		}
		colors.resize(parents.size());
		for (int i = 0; i < parents.size(); ++i) {
			colors[i] = compact_parents[get_parents(parents, i)];
		}
		printf("count %d\n", count);
		std::vector<std::vector<int> > neighbor_lists[4];
		for (int j = 0; j < 4; ++j)
			neighbor_lists[j].resize(num_v);
		Vector2i diffs[4] = { Vector2i(-1, 0), Vector2i(1, 0), Vector2i(0, -1), Vector2i(0, 1) };
		for (int i = 0; i < edge_to_constraints.size(); ++i) {
			if (cuts.count(edge_values[i]) == 0) {
				if (singularity_pos[get_parents(parents, edge_values[i].x)] || singularity_pos[get_parents(parents, edge_values[i].y)])
					continue;
				auto& edge_c = edge_to_constraints[i];
				int v0 = edge_c[0];
				int v1 = edge_c[2];

				int orientp0 = (get_parents_orient(parentss, v0) + edge_c[1]) % 4;
				int orientp1 = (get_parents_orient(parentss, v1) + edge_c[3]) % 4;
				int p1 = get_parents(parents, edge_values[i].x);
				int p2 = get_parents(parents, edge_values[i].y);
				if (p1 == p2)
					continue;

				auto diff = rshift90(edge_diff[i], orientp0);
				for (int j = 0; j < 4; ++j) {
					if (diffs[j] == diff)
						neighbor_lists[j][compact_parents[p1]].push_back(p2);
					if (diffs[j] == -diff)
						neighbor_lists[j][compact_parents[p2]].push_back(p1);
				}
			}
		}
		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < neighbor_lists[i].size(); ++j) {
				for (int k = 0; k < neighbor_lists[i][j].size(); ++k) {
					for (int l = k + 1; l < neighbor_lists[i][j].size(); ++l) {
						int p1 = get_parents(parents, neighbor_lists[i][j][k]);
						if (p1 != neighbor_lists[i][j][k])
							continue;
						int p2 = get_parents(parents, neighbor_lists[i][j][l]);
						if (p2 != neighbor_lists[i][j][l])
							continue;
						if (p1 == p2)
							continue;
						if (ranks[p1] > ranks[p2]) {
							ranks[p1] += ranks[p2];
							parents[p2] = std::make_pair(p1, 0);
							update = true;
						}
						else {
							ranks[p2] += ranks[p1];
							parents[p1] = std::make_pair(p2, 0);
							update = true;
						}
					}
				}
			}
		}
		if (!update)
			break;
	}


	/*
	int ttt = 0;
	for (int thres = 1; thres <= 0; ++thres) {
		bool update = true;
		while (update) {
			num_v = 0;
			for (int i = 0; i < parents.size(); ++i) {
				if (parents[i].first == i) {
					compact_parents[i] = num_v++;
				}
			}
			colors.resize(parents.size());
			for (int i = 0; i < parents.size(); ++i) {
				colors[i] = compact_parents[get_parents(parents, i)];
			}


			update = false;
			std::vector<std::map<int, int> > degrees(num_v);
			for (int i = 0; i < edge_values.size(); ++i) {
				int v1 = colors[edge_values[i].x];
				int v2 = colors[edge_values[i].y];
				if (v1 != v2) {
					degrees[v1][v2] = i;
					degrees[v2][v1] = i;
				}
			}
			
			std::vector<int> visited(num_v, 0);
			for (int i = 0; i < edge_diff.size(); ++i) {
				int p1 = colors[edge_values[i].x];
				int p2 = colors[edge_values[i].y];
				if (visited[p1] || visited[p2]) {
					continue;
				}
				if ((degrees[p1].size() <= thres || degrees[p2].size() <= thres) && p1 != p2) {
					if (thres == 2) {
						printf("%d %d\n", p1, p2);
						ttt += 1;
					}
					visited[p1] = 1;
					visited[p2] = 1;
					int vv0 = edge_values[i].x;
					int vv1 = edge_values[i].y;
					int v0 = get_parents(parents, edge_values[i].x);
					int v1 = get_parents(parents, edge_values[i].y);
					auto index = compat_orientation_extrinsic_index_4(Q.col(edge_values[i].x),
						N.col(edge_values[i].x), Q.col(edge_values[i].y), N.col(edge_values[i].y));
					int diff = (index.second - index.first + 4) % 4;

					if (v0 == v1) {
						int t = get_parents_orient(parents, vv0) - get_parents_orient(parents, vv1);
						if (diff != (4 + get_parents_orient(parents, vv0) - get_parents_orient(parents, vv1)) % 4) {
						}
						continue;
					}
					if (degrees[p2].size() <= thres) {
						ranks[v0] += ranks[v1];
						diff = (8 - diff - get_parents_orient(parents, edge_values[i].y) + get_parents_orient(parents, edge_values[i].x)) % 4;
						parents[v1] = std::make_pair(v0, diff);
						CollapseVertex(edge_values[i].y, edge_values[i].x);
//						if (thres == 1)
							update = true;
///						else
//							break;
					}
					else {
						ranks[v1] += ranks[v0];
						diff = (4 - get_parents_orient(parents, edge_values[i].x) + diff + get_parents_orient(parents, edge_values[i].y)) % 4;
						parents[v0] = std::make_pair(v1, diff);
						CollapseVertex(edge_values[i].x, edge_values[i].y);
//						if (thres == 1)
							update = true;
//						else
//							break;
					}

				}
			}
		}
	}
//	CheckSanity(0);
	*/

	O_compact.resize(num_v, Vector3d::Zero());
	counter.resize(num_v, 0);
	std::vector<Vector3d> Q_compact(num_v, Vector3d::Zero());
	std::vector<Vector3d> N_compact(num_v, Vector3d::Zero());
	for (int i = 0; i < O.cols(); ++i) {
		O_compact[colors[i]] += O.col(i);
		N_compact[colors[i]] += N.col(i);
		Q_compact[colors[i]] += rotate90_by(Q.col(i), N.col(i), get_parents_orient(parents, i));
		counter[colors[i]] += 1;
	}
	for (int i = 0; i < O_compact.size(); ++i) {
		O_compact[i] /= counter[i];
		N_compact[i].normalize();
		Q_compact[i].normalize();
	}
	Es.clear();
	for (auto& p : singularities) {
		int v1 = colors[F(0, p.first)];
		int v2 = colors[F(1, p.first)];
		int v3 = colors[F(2, p.first)];
		int t = (v1 == v2) + (v1 == v3) + (v2 == v3);

		singular_patches.insert(v1);
		singular_patches.insert(v2);
		singular_patches.insert(v3);
		if (t == 1) {
			if (v1 != v2)
				singular_e.insert(DEdge(v1, v2));
			else if (v1 != v3)
				singular_e.insert(DEdge(v1, v3));
			else if (v2 != v3)
				singular_e.insert(DEdge(v2, v3));
			else
				printf("!!!\n");
		}
		else {
			singular_e.insert(DEdge(v1, v2));
			singular_e.insert(DEdge(v2, v3));
			singular_e.insert(DEdge(v1, v3));
		}
	}


	for (int i = 0; i < edge_diff.size(); ++i) {
		int v1 = edge_values[i].x;
		int v2 = edge_values[i].y;
		int p1 = colors[v1];
		int p2 = colors[v2];
		if (p1 == p2)
			continue;
		auto index = compat_orientation_extrinsic_index_4(Q.col(v1), N.col(v1), Q.col(v2), N.col(v2));
		int orient_diff = (-get_parents_orient(parents, v1) + get_parents_orient(parents, v2) + index.second - index.first + 8) % 4;
		if (singular_patches.count(p1))
			orient_diff = 0;
		Vector2i diffs[2];
		for (int j = 0; j < 2; ++j) {
			if (singular_patches.count(p1) == 0 && singular_e.count(DEdge(p1, p2)) == 0) {
				int orient = get_parents_orient(parents, v1);
				int orient_diff_current = (j == 0) ? orient_diff : (4 - orient_diff) % 4;
				auto diff = (j == 0) ? edge_diff[i] :
					rshift90(-edge_diff[i], (index.second - index.first + 4) % 4);
				diff = rshift90(diff, orient);
				diffs[j] = diff;
				if (Es.count(Edge(p1, p2)) != 0) {
					if (Es[Edge(p1, p2)].second != diff)
						Es[Edge(p1, p2)].second = Vector2i(10000, 10000);
					else if (Es[Edge(p1, p2)].first != orient_diff_current) {
						auto it = Es.find(Edge(p2, p1));
						if (it != Es.end())
							it->second.second = Vector2i(10000, 10000);
					}
				}
				else
					Es[Edge(p1, p2)] = std::make_pair(orient_diff_current, diff);
			}
			std::swap(v1, v2);
			std::swap(p1, p2);
		}
	}
	
	while (!Es.empty() && Es.begin()->second.second == Vector2i(10000, 10000)) {
		singular_e.insert(DEdge(Es.begin()->first.first, Es.begin()->first.second));
		Es.erase(Es.begin());
	}
	for (auto it = Es.begin(), prev_it = it; it != Es.end(); prev_it = it, ++it) {
		if (it->second.second == Vector2i(10000, 10000)) {
			singular_e.insert(DEdge(it->first.first, it->first.second));
			Es.erase(it);
			it = prev_it;
		}
	}
	
	/*
	for (auto& e : Es) {
		if (singular_patches.count(e.first.first) == 0 && singular_patches.count(e.first.second)) {
			printf("found %d %d\n", e.first.first, e.first.second);
			system("pause");
		}
	}
	for (auto it = Es.begin(); it != Es.end(); ++it) {
		if (abs(it->second.second[0]) + abs(it->second.second[1]) > 1) {
			int directions = Es.count(Edge(it->first.second, it->first.first));
			if (it->first.first < it->first.second || directions == 0) {
				int prev_v = it->first.first;
				int dir_x = it->second.second[0] > 0 ? 1 : -1;
				int dir_y = it->second.second[1] > 0 ? 1 : -1;
				Vector3d p = O_compact[it->first.first];
				Vector3d d = O_compact[it->first.second] - p;
				Vector3d& q = Q_compact[it->first.first];
				Vector3d& n = N_compact[it->first.first];
				Vector3d q_y = n.cross(q);
				Vector3d tx = (d.dot(q) / abs(it->second.second[0])) * q;
				Vector3d ty = (d.dot(q_y) / abs(it->second.second[1])) * q_y;
				Vector3d tn = (d.dot(n) / (abs(it->second.second[0]) + abs(it->second.second[1]))) * n;
				double weight = counter[it->first.first];
				double dweight = (counter[it->first.second] - counter[it->first.first]) / (abs(it->second.second[0]) + abs(it->second.second[1]));
				for (int i = 0; i != it->second.second[0]; i += dir_x) {
					int next_v;
					int diff = 0;
					int new_v = 0;
					Vector2i reverse_dir = Vector2i(-dir_x, 0);
					p += tx + tn;
					weight += dweight;
					if (it->second.second[1] == 0 && i + dir_x == it->second.second[0]) {
						next_v = it->first.second;
						diff = it->second.first;
						reverse_dir = rshift90(reverse_dir, diff);
					}
					else {
						next_v = num_v++;
						if (next_v == 22419)
							next_v = next_v;
						O_compact.push_back(p);
						counter.push_back(dweight);
						new_v = 1;
					}
					Es[Edge(prev_v, next_v)] = std::make_pair(diff, Vector2i(dir_x, 0));
					if (directions || new_v) {
						Es[Edge(next_v, prev_v)] = std::make_pair((4 - diff) % 4, reverse_dir);
					}
					prev_v = next_v;
				}
				for (int i = 0; i != it->second.second[1]; i += dir_y) {
					int next_v;
					int diff = 0;
					int new_v = 0;
					Vector2i reverse_dir = Vector2i(0, -dir_y);
					p += ty + tn;
					weight += dweight;
					if (i + dir_y == it->second.second[0]) {
						next_v = it->first.second;
						diff = it->second.first;
						reverse_dir = rshift90(reverse_dir, diff);
					}
					else {
						next_v = num_v++;
						if (next_v == 22419)
							next_v = next_v;
						O_compact.push_back(p);
						counter.push_back(dweight);
						new_v = 1;
					}
					Es[Edge(prev_v, next_v)] = std::make_pair(diff, Vector2i(0, dir_y));
					if (directions || new_v) {
						Es[Edge(next_v, prev_v)] = std::make_pair((4 - diff) % 4, reverse_dir);
					}
					prev_v = next_v;
				}
			}
		}
	}*/
	

	/*
	while (!Es.empty() && abs(Es.begin()->second.second[0]) + abs(Es.begin()->second.second[1]) > 1)
		Es.erase(Es.begin());
	for (auto it = Es.begin(), prev_it = it; it != Es.end(); prev_it = it, ++it) {
		if (abs(it->second.second[0]) + abs(it->second.second[1]) > 1) {
			Es.erase(it);
			it = prev_it;
		}
	}
	*/
	shrink_parents.resize(num_v);
	shrink_compact_indices.resize(num_v);
	shrink_ranks.resize(num_v, 1);
	for (int i = 0; i < shrink_parents.size(); ++i) {
		shrink_parents[i] = std::make_pair(i, 0);
		shrink_compact_indices[i] = i;
	}
	printf("shrink_parents %d\n", shrink_parents.size());
	compact_num_v = num_v;
	int t = 0;
	
	while (true) {
		printf("update...\n");
		break;
		bool update = false;
		std::vector<std::map<int, std::pair<int, Vector2i> > > links(compact_num_v);
		for (auto& s : singular_patches)
			if (singular_patches.count(get_parents(shrink_parents, s)) == 0) {
				printf("fail %d %d!\n", s, get_parents(shrink_parents, s));
				system("pause");
			}
		for (auto& e : Es) {
			int v1 = shrink_compact_indices[get_parents(shrink_parents, e.first.first)];
			int v2 = get_parents(shrink_parents, e.first.second);
			if (singular_patches.count(v1) || singular_patches.count(v2))
				continue;
			int orient = get_parents_orient(shrink_parents, e.first.first);
			auto diff = rshift90(e.second.second, orient);
			// vi - vo = vi - e.first
			int diff_orient = (-orient + e.second.first + get_parents_orient(shrink_parents, e.first.second) + 4) % 4;
			if (links[v1].count(v2)) {
				if (links[v1][v2].second != diff || links[v1][v2].first != diff_orient)
					links[v1][v2].second = Vector2i(10000, 10000);
			}
			else {
				links[v1][v2] = std::make_pair(diff_orient, diff);
			}
		}
		Vector2i entries[] = { Vector2i(1, 0), Vector2i(0, 1), Vector2i(-1, 0), Vector2i(0, -1) };
		auto shrink_parents_buf = shrink_parents;
		std::set<DEdge> different_pairs;
		for (int i = 0; i < links.size(); ++i) {
			std::list<int> l[4];
			for (int j = 0; j < 4; ++j) {
				for (auto& e : links[i]) {
					if (e.second.second == entries[j]) {
						l[j].push_back(get_parents(shrink_parents, e.first));
					}
				}
			}
			for (int j = 0; j < 4; ++j) {
				for (auto& p1 : l[j]) {
					for (int k = j + 1; k < 4; ++k) {
						for (auto& p2 : l[k]) {
							different_pairs.insert(DEdge(p1, p2));
						}
					}
					different_pairs.insert(DEdge(p1, i));
				}
			}
		}
		for (int j = 0; j < 4; ++j) {
			for (int i = 0; i < links.size(); ++i) {
				std::list<std::pair<int, int> > l;
				for (auto& e : links[i]) {
					if (e.second.second == entries[j]) {
						l.push_back(std::make_pair(e.first, e.second.first));
					}
				}
				if (l.empty())
					continue;
				auto& v = l.front();
				for (auto& e : l) {
					int v0 = get_parents(shrink_parents, v.first);
					int v1 = get_parents(shrink_parents, e.first);
					if (v0 != v1) {
						if (singular_patches.count(v1) == 0 && (shrink_ranks[v0] > shrink_ranks[v1] || singular_patches.count(v0) != 0)) {
							shrink_ranks[v0] += shrink_ranks[v1];
							int diff = (get_parents_orient(shrink_parents, v0) + v.second - e.second - get_parents_orient(shrink_parents, v1) + 8) % 4;
							shrink_parents[v1] = std::make_pair(v0, diff);
							update = true;
						}
						else {
							shrink_ranks[v1] += shrink_ranks[v0];
							int diff = (get_parents_orient(shrink_parents, v1) + e.second - v.second - get_parents_orient(shrink_parents, v0) + 8) % 4;
							shrink_parents[v0] = std::make_pair(v1, diff);
							update = true;
						}
					}
					else {
						if (v.second != e.second) {
							printf("singular %d\n", v0);
							singular_patches.insert(v0);
						}
					}
				}
			}
		}
		if (!update)
			break;
		else {
			compact_num_v = 0;
			for (int i = 0; i < shrink_parents.size(); ++i) {
				if (shrink_parents[i].first == i)
					shrink_compact_indices[i] = compact_num_v++;
			}
		}
		t += 1;
	}

	for (auto& f : singularities) {
		printf("sing: %d %d %d %d\n", f.first,
			shrink_compact_indices[get_parents(shrink_parents, colors[F(0, f.first)])],
			shrink_compact_indices[get_parents(shrink_parents, colors[F(1, f.first)])],
			shrink_compact_indices[get_parents(shrink_parents, colors[F(2, f.first)])]);
	}

	UpdateMesh();
	CheckSanity(0);
	return;
	mE_extracted.clear();
	mO_extracted.resize(compact_num_v, Vector3d(0, 0, 0));
	std::vector<double> shrink_counter(compact_num_v, 0);
	for (int i = 0; i < O_compact.size(); ++i) {
		int p = shrink_compact_indices[get_parents(shrink_parents, i)];
		shrink_counter[p] += counter[i];
		mO_extracted[p] += O_compact[i] * counter[i];
	}
	for (int i = 0; i < mO_extracted.size(); ++i) {
		mO_extracted[i] /= shrink_counter[i];
	}
	std::set<DEdge> edges;
	for (auto& e : Es) {
		int v1 = shrink_compact_indices[get_parents(shrink_parents, e.first.first)];
		int v2 = shrink_compact_indices[get_parents(shrink_parents, e.first.second)];
		edges.insert(DEdge(v1, v2));
	}

	std::set<int> singular_patches_buf;
	std::set<DEdge> singular_e_buf;
	for (auto& s : singular_patches) {
		singular_patches_buf.insert(shrink_compact_indices[get_parents(shrink_parents, s)]);
	}
	for (auto& e : singular_e) {
		auto diff = edge_diff[edge_ids[e]];
		int v1 = shrink_compact_indices[get_parents(shrink_parents, e.x)];
		int v2 = shrink_compact_indices[get_parents(shrink_parents, e.y)];
		singular_e_buf.insert(DEdge(v1, v2));
	}
	std::swap(singular_e, singular_e_buf);
	std::swap(singular_patches, singular_patches_buf);
	mE_extracted.insert(mE_extracted.end(), edges.begin(), edges.end());
	printf("%d %d\n", mO_extracted.size(), mE_extracted.size());

}

void Parametrizer::BuildIntegerConstraints()
{
	auto& F = hierarchy.mF;
	auto& Q = hierarchy.mQ[0];
	auto& N = hierarchy.mN[0];
	std::vector<Vector2i> sign_indices;
	edge_to_constraints.resize(edge_ids.size(), Vector4i(-1,-1,-1,-1));
	for (int i = 0; i < F.cols(); ++i) {
		int v0 = F(0, i);
		int v1 = F(1, i);
		int v2 = F(2, i);
		DEdge e0(v0, v1), e1(v1, v2), e2(v2, v0);
		int eid[3] = { edge_ids[e0], edge_ids[e1], edge_ids[e2] };
		Vector2i vid[3];
		for (int i = 0; i < 3; ++i) {
			vid[i] = Vector2i(eid[i] * 2 + 1, eid[i] * 2 + 2);
		}
		auto index1 = compat_orientation_extrinsic_index_4(Q.col(v0), N.col(v0), Q.col(v1), N.col(v1));
		auto index2 = compat_orientation_extrinsic_index_4(Q.col(v0), N.col(v0), Q.col(v2), N.col(v2));
		int rank1 = (index1.first - index1.second + 4) % 4;
		int rank2 = (index2.first - index2.second + 4) % 4;
		int orients[3] = { 0 };
		if (v1 < v0) {
			vid[0] = -rshift90(vid[0], rank1);
			orients[0] = (rank1 + 2) % 4;
		}
		if (v2 < v1) {
			vid[1] = -rshift90(vid[1], rank2);
			orients[1] = (rank2 + 2) % 4;
		}
		else {
			vid[1] = rshift90(vid[1], rank1);
			orients[1] = rank1;
		}
		if (v2 < v0) {
			vid[2] = rshift90(vid[2], rank2);
			orients[2] = rank2;
		}
		else {
			vid[2] = -vid[2];
			orients[2] = 2;
		}

		edge_to_constraints[eid[0]][(v0 > v1) * 2] = i;
		edge_to_constraints[eid[0]][(v0 > v1) * 2 + 1] = orients[0];
		edge_to_constraints[eid[1]][(v1 > v2) * 2] = i;
		edge_to_constraints[eid[1]][(v1 > v2) * 2 + 1] = orients[1];
		edge_to_constraints[eid[2]][(v2 > v0) * 2] = i;
		edge_to_constraints[eid[2]][(v2 > v0) * 2 + 1] = orients[2];

		for (int k = 0; k < 3; ++k) {
			sign_indices.push_back(vid[k]);
		}
	}
	parentss.resize(F.cols());
	std::vector<int> ranks(F.cols(), 1);
	for (int i = 0; i < parentss.size(); ++i)
		parentss[i] = std::make_pair(i, 0);
	for (int i = 0; i < edge_to_constraints.size(); ++i) {
		auto& edge_c = edge_to_constraints[i];
		int v0 = edge_c[0];
		int v1 = edge_c[2];
		int orient1 = edge_c[1];
		int orient0 = (edge_c[3] + 2) % 4;

		int p0 = get_parents(parentss, v0);
		int p1 = get_parents(parentss, v1);
		int orientp0 = get_parents_orient(parentss, v0);
		int orientp1 = get_parents_orient(parentss, v1);

		if (p0 == p1) {
			continue;
		}
		if (ranks[p1] < ranks[p0]) {
			ranks[p0] += ranks[p1];
			parentss[p1].first = p0;
			parentss[p1].second = (orient1 - orient0 + orientp0 - orientp1 + 8) % 4;
		}
		else {
			ranks[p1] += ranks[p0];
			parentss[p0].first = p1;
			parentss[p0].second = (orient0 - orient1 + orientp1 - orientp0 + 8) % 4;
		}
	}

	std::vector<Vector3i> sing_diff;
	std::vector<std::map<int, std::pair<int, int> > > sing_maps;
	for (int i = 0; i < sign_indices.size(); i += 3) {
		int f = i / 3;
		int orient = get_parents_orient(parentss, f);		
		for (int j = 0; j < 3; ++j) {
			int v1 = F(j, f);
			int v2 = F((j + 1) % 3, f);
			int eid = edge_ids[DEdge(v1, v2)];
			sign_indices[i + j] = rshift90(sign_indices[i + j], orient);
		}
		for (int j = 0; j < 2; ++j) {
			Vector3i sign, ind;
			for (int k = 0; k < 3; ++k) {
				ind[k] = abs(sign_indices[i + k][j]);
				sign[k] = sign_indices[i + k][j] / ind[k];
				ind[k] -= 1;
			}
			constraints_index.push_back(ind);
			constraints_sign.push_back(sign);
		}
		if (singularities.count(f)) {
			int orient_base = singularities[f];
			Vector3i diffs;
			for (int j = 0; j < 3; ++j) {
				auto sign_index = rshift90(sign_indices[i + j], orient_base);
				int total_diff = 0;
				for (int k = 0; k < 2; ++k) {
					auto ind = constraints_index[f * 2 + k];
					auto sign = constraints_sign[f * 2 + k];
					ind[j] = abs(sign_index[k]);
					sign[j] = sign_index[k] / abs(sign_index[k]);
					ind[j] -= 1;
					int diff1 = edge_diff[ind[0] / 2][ind[0] % 2];
					int diff2 = edge_diff[ind[1] / 2][ind[1] % 2];
					int diff3 = edge_diff[ind[2] / 2][ind[2] % 2];
					int diff = sign[0] * diff1 + sign[1] * diff2 + sign[2] * diff3;
					total_diff += diff;
				}
				diffs[j] = total_diff;
			}
			sing_diff.push_back(diffs);
		}
	}

	int total_flow = 0;
	for (int i = 0; i < constraints_index.size(); ++i) {
		if (singularities.count(i / 2)) {
			continue;
		}
		auto index = constraints_index[i];
		auto sign = constraints_sign[i];
		int diff1 = edge_diff[index[0] / 2][index[0] % 2];
		int diff2 = edge_diff[index[1] / 2][index[1] % 2];
		int diff3 = edge_diff[index[2] / 2][index[2] % 2];
		int diff = sign[0] * diff1 + sign[1] * diff2 + sign[2] * diff3;
		total_flow += diff;
	}
	sing_maps.resize(sing_diff.size() + 1);
	sing_maps[0][total_flow] = std::make_pair(0, 0);
	for (int i = 0; i < sing_diff.size(); ++i) {
		auto& prev = sing_maps[i];
		auto& next = sing_maps[i + 1];
		for (auto& p : prev) {
			for (int j = 0; j < 3; ++j) {
				int v = p.first + sing_diff[i][j];
				int t = p.second.first + abs(sing_diff[i][j]);
				auto it = next.find(v);
				if (it == next.end())
					next[v] = std::make_pair(t, j);
				else
					if (t < it->second.first)
						it->second = std::make_pair(t, j);
			}
		}
	}
	if (sing_maps.back().count(0) == 0) {
		printf("No Zero map!\n");
		exit(0);
	}
	std::vector<int> sing_selection;
	int target_flow = 0;
	for (int i = sing_diff.size(); i > 0; i--) {
		auto p = sing_maps[i][target_flow];
		target_flow -= sing_diff[i - 1][p.second];
		sing_selection.push_back(p.second);
	}
	std::reverse(sing_selection.begin(), sing_selection.end());
	int sing_count = 0;
	for (auto& f : singularities) {
		int select = sing_selection[sing_count++];
		auto& index1 = constraints_index[f.first * 2];
		auto& index2 = constraints_index[f.first * 2 + 1];
		auto& sign1 = constraints_sign[f.first * 2];
		auto& sign2 = constraints_sign[f.first * 2 + 1];
		auto diff = Vector2i(sign1[select] * (index1[select] + 1),
			sign2[select] * (index2[select] + 1));
		diff = rshift90(diff, f.second);
		index1[select] = abs(diff[0]);
		sign1[select] = diff[0] / abs(diff[0]);
		index1[select] -= 1;
		index2[select] = abs(diff[1]);
		sign2[select] = diff[1] / abs(diff[1]);
		index2[select] -= 1;
	}
	total_flow = 0;

	variables.resize(edge_diff.size() * 2, std::make_pair(Vector2i(-1, -1), 0));
	for (int i = 0; i < constraints_index.size(); ++i) {
		auto index = constraints_index[i];
		auto sign = constraints_sign[i];
		int diff1 = edge_diff[index[0] / 2][index[0] % 2];
		int diff2 = edge_diff[index[1] / 2][index[1] % 2];
		int diff3 = edge_diff[index[2] / 2][index[2] % 2];
		int diff = sign[0] * diff1 + sign[1] * diff2 + sign[2] * diff3;
		total_flow += diff;
		if (diff != 0 && pos_sing.count(i / 2) == 0 && singularities.count(i / 2) == 0) {
			printf("FAIL!\n");
			system("pause");
		}
		for (int j = 0; j < 3; ++j) {
			auto& p = variables[index[j]].first;
			if (sign[j] > 0)
				p[0] = i;
			else
				p[1] = i;
			variables[index[j]].second += sign[j];
		}
	}
	cuts.clear();
	for (int i = 0; i < variables.size(); ++i) {
		if (variables[i].second != 0) {
			cuts.insert(edge_values[i / 2]);
		}
	}

	ComputeMaxFlow();
}

void Parametrizer::ComputeMaxFlow()
{
	//"p min 6 8
	//c min-cost flow problem with 6 nodes and 8 arcs
	//n 1 10
	//c supply of 10 at node 1
	//n 6 -10
	//c demand of 10 at node 6
	//c arc list follows
	//c arc has <tail> <head> <capacity l.b.> <capacity u.b> <cost>
	//a 1 2 0 4 1
	//a 1 3 0 8 5
	//a 2 3 0 5 0
	//a 3 5 0 10 1
	//a 5 4 0 8 0
	//a 5 6 0 8 9
	//a 4 2 0 8 1
	//a 4 6 0 8 1"
	auto& Q = hierarchy.mQ[0];
	auto& N = hierarchy.mN[0];
	auto& V = hierarchy.mV[0];
	int num_nodes = constraints_index.size() + 2;
	std::vector<std::pair<Vector2i, int> > arcs;
	std::vector<std::pair<int, int> > arc_ids;
	for (int i = 0; i < variables.size(); ++i) {
		if (variables[i].second == 0) {
			int current_v = edge_diff[i / 2][i % 2];
			int v1 = edge_values[i / 2].x;
			int v2 = edge_values[i / 2].y;
			arcs.push_back(std::make_pair(variables[i].first, current_v));
			Vector3d q = Q.col(v1);
			Vector3d d = V.col(v2) - V.col(v1);
			Vector3d n = N.col(v1);
			double t = 0;
			if (i % 2 == 0)
				t = q.dot(d);
			else
				t = d.dot(n.cross(d));
			arc_ids.push_back(std::make_pair(i, t > 0 ? 1 : -1));
		}
	}
	int supply = 0;
	for (int i = 0; i < constraints_index.size(); ++i) {
		int diff = 0;
		for (int j = 0; j < 3; ++j) {
			int ind = constraints_index[i][j];
			diff += constraints_sign[i][j] * edge_diff[ind / 2][ind % 2];
		}
		if (diff > 0) {
			arcs.push_back(std::make_pair(Vector2i(-1, i), diff));
			supply += diff;
		}
		else {
			arcs.push_back(std::make_pair(Vector2i(i, constraints_index.size()), -diff));
		}
	}
	std::vector<std::map<int, std::pair<int, int> > > graph(constraints_index.size() + 2);
	std::map<Edge, std::pair<int, int> > edge_to_variable;
	for (int i = 0; i < arcs.size(); ++i) {
		int v1 = arcs[i].first[0] + 1;
		int v2 = arcs[i].first[1] + 1;
		int c = arcs[i].second;
		if (v1 == 0 || v2 == constraints_index.size() + 1) {
			graph[v1][v2] = std::make_pair(c, 0);
			graph[v2][v1] = std::make_pair(0, 0);
		}
		else {
			int t = arc_ids[i].second;
			graph[v1][v2] = std::make_pair(std::max(0, c + (t == -1 ? 1 : 1)), 0);
			graph[v2][v1] = std::make_pair(std::max(0, -c + (t == 1 ? 1 : 1)), 0);
			edge_to_variable[Edge(v1, v2)] = std::make_pair(arc_ids[i].first, -1);
			edge_to_variable[Edge(v2, v1)] = std::make_pair(arc_ids[i].first, 1);
		}
	}

	int total_flow = 0;
	while (true) {
		std::vector<int> visited(graph.size(), 0);
		std::vector<std::pair<int, int> > q;
		q.push_back(std::make_pair(0, -1));
		visited[0] = 1;
		int front = 0;
		bool found = false;
		while (front < q.size() && !found) {
			for (auto& p : graph[q[front].first]) {
				if (visited[p.first] == 0 && p.second.first - p.second.second > 0) {
					visited[p.first] = 1;
					q.push_back(std::make_pair(p.first, front));
					if (p.first == constraints_index.size() + 1) {
						found = true;
						break;
					}
				}
			}
			front += 1;
		}
		if (found) {
			int flows = 0x7fffffff;
			int r = q.size() - 1;
			while (q[r].second != -1) {
				int v2 = q[r].first;
				int v1 = q[q[r].second].first;
				auto& p = graph[v1][v2];
				flows = std::min(flows, p.first - p.second);
				r = q[r].second;
			}
			r = q.size() - 1;
			while (q[r].second != -1) {
				int v2 = q[r].first;
				int v1 = q[q[r].second].first;
				graph[v1][v2].second += flows;
				graph[v2][v1].second -= flows;
				r = q[r].second;
			}
			total_flow += flows;
		}
		else
			break;
	}
	printf("supply and flow: %d %d\n", supply, total_flow);

	for (int i = 1; i < graph.size() - 1; ++i) {
		for (auto& p : graph[i]) {
			if (p.first >= 1 && p.first < graph.size() - 1) {
				int flow = p.second.second;
				if (flow > 0) {
					auto q = edge_to_variable[Edge(i, p.first)];
					edge_diff[q.first / 2][q.first % 2] += q.second * flow;
					if (std::abs(edge_diff[q.first / 2][q.first % 2]) > 1) {
						printf("FAIL!\n");
					}
				}
			}
		}
	}
}