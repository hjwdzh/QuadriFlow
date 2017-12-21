#include "Parametrizer.h"

#include "loader.h"
#include "MergeVertex.h"
#include "subdivide.h"
#include "dedge.h"
#include "AdjacentMatrix.h"
#include "field_math.h"

#include <fstream>
#include <queue>

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
bool remove_unnecessary_edges) {

	double scale = mRes.mScale, inv_scale = 1 / scale;

	auto compat_orientation = compat_orientation_extrinsic_4;
	auto compat_position = compat_position_extrinsic_index_4;

	const MatrixXd &Q = mRes.mQ[0], &O = mRes.mO[0], &N = mRes.mN[0], &V = mRes.mV[0];
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

			double error = 0;
			std::pair<Vector2i, Vector2i> shift = compat_position(
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
	const MatrixXd &N = hierarchy.mN[0], &Q = hierarchy.mQ[0];
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

void Parametrizer::ExtractMesh() {
	printf("extract_graph\n");
	std::map<uint32_t, uint32_t> vertex_map;
	extract_graph(hierarchy, adj_extracted,
		mV_extracted, mN_extracted, vertex_map, true, true);
	MatrixXi F_extr;
	MatrixXd Nf_extr;

	auto adj_new = adj_extracted;
	extract_faces(adj_new, mV_extracted, mN_extracted, Nf_extr, F_extr,
		hierarchy.mScale, true);

	write_obj("result.obj", F_extr, mV_extracted, MatrixXd(), Nf_extr, MatrixXd(), MatrixXd());
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

void Parametrizer::SaveToFile(FILE* fp) {
	Save(fp, vertex_singularities);
	Save(fp, singularities);
	Save(fp, V);
	Save(fp, N);
	Save(fp, Nf);
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