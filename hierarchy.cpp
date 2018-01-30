#include "hierarchy.h"
#include "field_math.h"
#include <fstream>
#include "tbb_common.h"
#include "pss/parallel_stable_sort.h"
#include "pcg32/pcg32.h"
Hierarchy::Hierarchy()
{
	mAdj.resize(MAX_DEPTH + 1);
	mV.resize(MAX_DEPTH + 1);
	mN.resize(MAX_DEPTH + 1);
	mA.resize(MAX_DEPTH + 1);
	mPhases.resize(MAX_DEPTH + 1);
	mToLower.resize(MAX_DEPTH);
	mToUpper.resize(MAX_DEPTH);
}

void Hierarchy::Initialize(double scale, int with_scale)
{
	this->with_scale = with_scale;
	generate_graph_coloring_deterministic(mAdj[0], mV[0].cols(), mPhases[0]);
	for (int i = 0; i < MAX_DEPTH; ++i) {
		DownsampleGraph(mAdj[i], mV[i], mN[i], mA[i], mV[i + 1], mN[i + 1], mA[i + 1],
			mToUpper[i], mToLower[i], mAdj[i + 1]);
		generate_graph_coloring_deterministic(mAdj[i + 1], mV[i + 1].cols(), mPhases[i + 1]);
		if (mV[i + 1].cols() == 1) {
			mAdj.resize(i + 2);
			mV.resize(i + 2);
			mN.resize(i + 2);
			mA.resize(i + 2);
			mToUpper.resize(i + 1);
			mToLower.resize(i + 1);
			break;
		}
	}
	mQ.resize(mV.size());
	mO.resize(mV.size());

	if (with_scale) {
		mS.resize(mV.size());
		mK.resize(mV.size());
	}

	mScale = scale;
#pragma omp parallel for
	for (int i = 0; i<mV.size(); ++i) {
		mQ[i].resize(mN[i].rows(), mN[i].cols());
		mO[i].resize(mN[i].rows(), mN[i].cols());
		if (with_scale) {
			mS[i].resize(2, mN[i].cols());
			mK[i].resize(2, mN[i].cols());
		}
		for (int j = 0; j < mN[i].cols(); ++j) {
			Vector3d s, t;
			coordinate_system(mN[i].col(j), s, t);
			double angle = ((double)rand()) / RAND_MAX * 2 * M_PI;
			double x = ((double)rand()) / RAND_MAX * 2 - 1.f;
			double y = ((double)rand()) / RAND_MAX * 2 - 1.f;
			mQ[i].col(j) = s * std::cos(angle) + t * std::sin(angle);
			mO[i].col(j) = mV[i].col(j) + (s * x + t * y) * scale;
			if (with_scale) {
				mS[i].col(j) = Vector2d(1.0f, 1.0f);
				mK[i].col(j) = Vector2d(0.0, 0.0);
			}
		}
	}
#ifdef WITH_CUDA
	printf("copy to device...\n");
	CopyToDevice();
	printf("copy to device finish...\n");
#endif
}

void Hierarchy::generate_graph_coloring_deterministic(const AdjacentMatrix &adj, int size,
	std::vector<std::vector<int> > &phases) {
	struct ColorData {
		uint8_t nColors;
		uint32_t nNodes[256];
		ColorData() : nColors(0) { }
	};

	const uint8_t INVALID_COLOR = 0xFF;
	phases.clear();

	/* Generate a permutation */
	std::vector<uint32_t> perm(size);
	std::vector<tbb::spin_mutex> mutex(size);
	for (uint32_t i = 0; i<size; ++i)
		perm[i] = i;

	tbb::parallel_for(
		tbb::blocked_range<uint32_t>(0u, size, GRAIN_SIZE),
		[&](const tbb::blocked_range<uint32_t> &range) {
		pcg32 rng;
		rng.advance(range.begin());
		for (uint32_t i = range.begin(); i != range.end(); ++i) {
			uint32_t j = i, k =
				rng.nextUInt(size - i) + i;
			if (j == k)
				continue;
			if (j > k)
				std::swap(j, k);
			tbb::spin_mutex::scoped_lock l0(mutex[j]);
			tbb::spin_mutex::scoped_lock l1(mutex[k]);
			std::swap(perm[j], perm[k]);
		}
	}
	);

	std::vector<uint8_t> color(size, INVALID_COLOR);
	ColorData colorData = tbb::parallel_reduce(
		tbb::blocked_range<uint32_t>(0u, size, GRAIN_SIZE),
		ColorData(),
		[&](const tbb::blocked_range<uint32_t> &range, ColorData colorData) -> ColorData {
		std::vector<uint32_t> neighborhood;
		bool possible_colors[256];

		for (uint32_t pidx = range.begin(); pidx != range.end(); ++pidx) {
			uint32_t i = perm[pidx];

			neighborhood.clear();
			neighborhood.push_back(i);
//			for (const Link *link = adj[i]; link != adj[i + 1]; ++link)
			for (auto& link : adj[i])
				neighborhood.push_back(link.id);
			std::sort(neighborhood.begin(), neighborhood.end());
			for (uint32_t j : neighborhood)
				mutex[j].lock();

			std::fill(possible_colors, possible_colors + colorData.nColors, true);

//			for (const Link *link = adj[i]; link != adj[i + 1]; ++link) {
			for (auto& link : adj[i]) {
				uint8_t c = color[link.id];
				if (c != INVALID_COLOR) {
					while (c >= colorData.nColors) {
						possible_colors[colorData.nColors] = true;
						colorData.nNodes[colorData.nColors] = 0;
						colorData.nColors++;
					}
					possible_colors[c] = false;
				}
			}

			uint8_t chosen_color = INVALID_COLOR;
			for (uint8_t j = 0; j<colorData.nColors; ++j) {
				if (possible_colors[j]) {
					chosen_color = j;
					break;
				}
			}
			if (chosen_color == INVALID_COLOR) {
				if (colorData.nColors == INVALID_COLOR - 1)
					throw std::runtime_error("Ran out of colors during graph coloring! "
					"The input mesh is very likely corrupt.");
				colorData.nNodes[colorData.nColors] = 1;
				color[i] = colorData.nColors++;
			}
			else {
				colorData.nNodes[chosen_color]++;
				color[i] = chosen_color;
			}

			for (uint32_t j : neighborhood)
				mutex[j].unlock();
		}
		return colorData;
	},
		[](ColorData c1, ColorData c2) -> ColorData {
		ColorData result;
		result.nColors = max(c1.nColors, c2.nColors);
		memset(result.nNodes, 0, sizeof(uint32_t) * result.nColors);
		for (uint8_t i = 0; i<c1.nColors; ++i)
			result.nNodes[i] += c1.nNodes[i];
		for (uint8_t i = 0; i<c2.nColors; ++i)
			result.nNodes[i] += c2.nNodes[i];
		return result;
	}
	);

	phases.resize(colorData.nColors);
	for (int i = 0; i<colorData.nColors; ++i)
		phases[i].reserve(colorData.nNodes[i]);

	for (uint32_t i = 0; i<size; ++i)
		phases[color[i]].push_back(i);

}


void Hierarchy::DownsampleGraph(const AdjacentMatrix adj, const MatrixXd &V,
	const MatrixXd &N, const VectorXd &A,
	MatrixXd &V_p, MatrixXd &N_p, VectorXd &A_p,
	MatrixXi &to_upper, VectorXi &to_lower,
	AdjacentMatrix& adj_p) {
	struct Entry {
		int i, j;
		double order;
		inline Entry() { i = j = -1; };
		inline Entry(int i, int j, double order) : i(i), j(j), order(order) { }
		inline bool operator<(const Entry &e) const { return order > e.order; }
		inline bool operator==(const Entry &e) const { return order == e.order; }
	};

	int nLinks = 0;
	for (auto& adj_i : adj)
		nLinks += adj_i.size();
	std::vector<Entry> entries(nLinks);
	std::vector<int> bases(adj.size());
	for (int i = 1; i < bases.size(); ++i) {
		bases[i] = bases[i - 1] + adj[i - 1].size();
	}

#pragma omp parallel for
	for (int i = 0; i < V.cols(); ++i) {
		int num = adj[i].size();
		int base = bases[i];
		auto& ad = adj[i];
		auto entry_it = entries.begin() + base;
		for (auto it = ad.begin(); it != ad.end(); ++it, ++entry_it) {
			int k = it->id;
			double dp = N.col(i).dot(N.col(k));
			double ratio = A[i] > A[k] ? (A[i] / A[k]) : (A[k] / A[i]);
			*entry_it = Entry(i, k, dp * ratio);
		}
	}

	pss::parallel_stable_sort(entries.begin(), entries.end(), std::less<Entry>());

	std::vector<bool> mergeFlag(V.cols(), false);

	int nCollapsed = 0;
	for (int i = 0; i<nLinks; ++i) {
		const Entry &e = entries[i];
		if (mergeFlag[e.i] || mergeFlag[e.j])
			continue;
		mergeFlag[e.i] = mergeFlag[e.j] = true;
		entries[nCollapsed++] = entries[i];
	}
	int vertexCount = V.cols() - nCollapsed;

	// Allocate memory for coarsened graph 
	V_p.resize(3, vertexCount);
	N_p.resize(3, vertexCount);
	A_p.resize(vertexCount);
	to_upper.resize(2, vertexCount);
	to_lower.resize(V.cols());

#pragma omp parallel for
	for (int i = 0; i < nCollapsed; ++i) {
		const Entry &e = entries[i];
		const double area1 = A[e.i], area2 = A[e.j], surfaceArea = area1 + area2;
		if (surfaceArea > RCPOVERFLOW)
			V_p.col(i) = (V.col(e.i) * area1 + V.col(e.j) * area2) / surfaceArea;
		else
			V_p.col(i) = (V.col(e.i) + V.col(e.j)) * 0.5f;
		Vector3d normal = N.col(e.i) * area1 + N.col(e.j) * area2;
		double norm = normal.norm();
		N_p.col(i) = norm > RCPOVERFLOW ? Vector3d(normal / norm)
			: Vector3d::UnitX();
		A_p[i] = surfaceArea;
		to_upper.col(i) << e.i, e.j;
		to_lower[e.i] = i; to_lower[e.j] = i;
	}

	int offset = nCollapsed;

	for (int i = 0; i < V.cols(); ++i) {
		if (!mergeFlag[i]) {
			int idx = offset++;
			V_p.col(idx) = V.col(i);
			N_p.col(idx) = N.col(i);
			A_p[idx] = A[i];
			to_upper.col(idx) << i, -1;
			to_lower[i] = idx;
		}
	}

	adj_p.resize(V_p.cols());
	std::vector<int> capacity(V_p.cols());
	std::vector<std::vector<Link> > scratches(V_p.cols());
#pragma omp parallel for
	for (int i = 0; i < V_p.cols(); ++i) {
		int t = 0;
		for (int j = 0; j < 2; ++j) {
			int upper = to_upper(j, i);
			if (upper == -1)
				continue;
			t += adj[upper].size();
		}
		scratches[i].reserve(t);
		adj_p[i].reserve(t);
	}
#pragma omp parallel for
	for (int i = 0; i < V_p.cols(); ++i) {
		auto& scratch = scratches[i];
		for (int j = 0; j<2; ++j) {
			int upper = to_upper(j, i);
			if (upper == -1)
				continue;
			auto& ad = adj[upper];
			for (auto& link : ad)
				scratch.push_back(Link(to_lower[link.id], link.weight));
		}
		std::sort(scratch.begin(), scratch.end());
		int id = -1;
		auto& ad = adj_p[i];
		for (auto& link : scratch) {
			if (link.id != i) {
				if (id != link.id) {
					ad.push_back(link);
					id = link.id;
				}
				else {
					ad.back().weight += link.weight;
				}
			}
		}
	}
}


void Hierarchy::SaveToFile(FILE* fp)
{
	Save(fp, mScale);
	Save(fp, mF);
	Save(fp, mE2E);
	Save(fp, mAdj);
	Save(fp, mV);
	Save(fp, mN);
	Save(fp, mA);
	Save(fp, mToLower);
	Save(fp, mToUpper);
	Save(fp, mQ);
	Save(fp, mO);
	Save(fp, mS);
	Save(fp, mK);
	Save(fp, this->mPhases);
}

void Hierarchy::LoadFromFile(FILE* fp)
{
	Read(fp, mScale);
	Read(fp, mF);
	Read(fp, mE2E);
	Read(fp, mAdj);
	Read(fp, mV);
	Read(fp, mN);
	Read(fp, mA);
	Read(fp, mToLower);
	Read(fp, mToUpper);
	Read(fp, mQ);
	Read(fp, mO);
	Read(fp, mS);
	Read(fp, mK);
	Read(fp, this->mPhases);
}

#ifdef WITH_CUDA
#include <cuda_runtime.h>

void Hierarchy::CopyToDevice()
{
	if (cudaAdj.empty()) {
		cudaAdj.resize(mAdj.size());
		cudaAdjOffset.resize(mAdj.size());
		for (int i = 0; i < mAdj.size(); ++i) {
			std::vector<int> offset(mAdj[i].size() + 1, 0);
			for (int j = 0; j < mAdj[i].size(); ++j) {
				offset[j + 1] = offset[j] + mAdj[i][j].size();
			}
			cudaMalloc(&cudaAdjOffset[i], sizeof(int) * (mAdj[i].size() + 1));
			cudaMemcpy(cudaAdjOffset[i], offset.data(), sizeof(int) * (mAdj[i].size() + 1), cudaMemcpyHostToDevice);
//			cudaAdjOffset[i] = (int*)malloc(sizeof(int) * (mAdj[i].size() + 1));
//			memcpy(cudaAdjOffset[i], offset.data(), sizeof(int) * (mAdj[i].size() + 1));

			cudaMalloc(&cudaAdj[i], sizeof(Link) * offset.back());
//			cudaAdj[i] = (Link*)malloc(sizeof(Link) * offset.back());
			std::vector<Link> plainlink(offset.back());
			for (int j = 0; j < mAdj[i].size(); ++j) {
				memcpy(plainlink.data() + offset[j], mAdj[i][j].data(), mAdj[i][j].size() * sizeof(Link));
			}
			cudaMemcpy(cudaAdj[i], plainlink.data(), plainlink.size() * sizeof(Link), cudaMemcpyHostToDevice);
		}
	}

	if (cudaN.empty()) {
		cudaN.resize(mN.size());
		for (int i = 0; i < mN.size(); ++i) {
			cudaMalloc(&cudaN[i], sizeof(glm::dvec3) * mN[i].cols());
//			cudaN[i] = (glm::dvec3*)malloc(sizeof(glm::dvec3) * mN[i].cols());
		}
	}
	for (int i = 0; i < mN.size(); ++i) {
		cudaMemcpy(cudaN[i], mN[i].data(), sizeof(glm::dvec3) * mN[i].cols(), cudaMemcpyHostToDevice);
//		memcpy(cudaN[i], mN[i].data(), sizeof(glm::dvec3) * mN[i].cols());
	}

	if (cudaV.empty()) {
		cudaV.resize(mV.size());
		for (int i = 0; i < mV.size(); ++i) {
			cudaMalloc(&cudaV[i], sizeof(glm::dvec3) * mV[i].cols());
//			cudaV[i] = (glm::dvec3*)malloc(sizeof(glm::dvec3) * mV[i].cols());
		}
	}
	for (int i = 0; i < mV.size(); ++i) {
		cudaMemcpy(cudaV[i], mV[i].data(), sizeof(glm::dvec3) * mV[i].cols(), cudaMemcpyHostToDevice);
//		memcpy(cudaV[i], mV[i].data(), sizeof(glm::dvec3) * mV[i].cols());
	}

	if (cudaQ.empty()) {
		cudaQ.resize(mQ.size());
		for (int i = 0; i < mQ.size(); ++i) {
			cudaMalloc(&cudaQ[i], sizeof(glm::dvec3) * mQ[i].cols());
//			cudaQ[i] = (glm::dvec3*)malloc(sizeof(glm::dvec3) * mQ[i].cols());
		}
	}
	for (int i = 0; i < mQ.size(); ++i) {
		cudaMemcpy(cudaQ[i], mQ[i].data(), sizeof(glm::dvec3) * mQ[i].cols(), cudaMemcpyHostToDevice);
//		memcpy(cudaQ[i], mQ[i].data(), sizeof(glm::dvec3) * mQ[i].cols());
	}
	if (cudaO.empty()) {
		cudaO.resize(mO.size());
		for (int i = 0; i < mO.size(); ++i) {
			cudaMalloc(&cudaO[i], sizeof(glm::dvec3) * mO[i].cols());
//			cudaO[i] = (glm::dvec3*)malloc(sizeof(glm::dvec3) * mO[i].cols());
		}
	}
	for (int i = 0; i < mO.size(); ++i) {
		cudaMemcpy(cudaO[i], mO[i].data(), sizeof(glm::dvec3) * mO[i].cols(), cudaMemcpyHostToDevice);
//		memcpy(cudaO[i], mO[i].data(), sizeof(glm::dvec3) * mO[i].cols());
	}
	if (cudaPhases.empty()) {
		cudaPhases.resize(mPhases.size());
		for (int i = 0; i < mPhases.size(); ++i) {
			cudaPhases[i].resize(mPhases[i].size());
			for (int j = 0; j < mPhases[i].size(); ++j) {
				cudaMalloc(&cudaPhases[i][j], sizeof(int) * mPhases[i][j].size());
//				cudaPhases[i][j] = (int*)malloc(sizeof(int) * mPhases[i][j].size());
			}
		}
	}
	for (int i = 0; i < mPhases.size(); ++i) {
		for (int j = 0; j < mPhases[i].size(); ++j) {
			cudaMemcpy(cudaPhases[i][j], mPhases[i][j].data(), sizeof(int) * mPhases[i][j].size(), cudaMemcpyHostToDevice);
//			memcpy(cudaPhases[i][j], mPhases[i][j].data(), sizeof(int) * mPhases[i][j].size());
		}
	}
	if (cudaToUpper.empty()) {
		cudaToUpper.resize(mToUpper.size());
		for (int i = 0; i < mToUpper.size(); ++i) {
			cudaMalloc(&cudaToUpper[i], mToUpper[i].cols() * sizeof(glm::ivec2));
//			cudaToUpper[i] = (glm::ivec2*)malloc(mToUpper[i].cols() * sizeof(glm::ivec2));
		}
	}
	for (int i = 0; i < mToUpper.size(); ++i) {
		cudaMemcpy(cudaToUpper[i], mToUpper[i].data(), sizeof(glm::ivec2) * mToUpper[i].cols(), cudaMemcpyHostToDevice);
//		memcpy(cudaToUpper[i], mToUpper[i].data(), sizeof(glm::ivec2) * mToUpper[i].cols());
	}
	cudaDeviceSynchronize();
}


void Hierarchy::CopyToHost()
{

}


#endif
