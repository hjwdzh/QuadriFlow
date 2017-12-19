#include "hierarchy.h"
#include "field_math.h"
#include <fstream>

Hierarchy::Hierarchy()
{
	mAdj.resize(MAX_DEPTH + 1);
	mV.resize(MAX_DEPTH + 1);
	mN.resize(MAX_DEPTH + 1);
	mA.resize(MAX_DEPTH + 1);
	mToLower.resize(MAX_DEPTH);
	mToUpper.resize(MAX_DEPTH);
}

void Hierarchy::Initialize(float scale)
{
	for (int i = 0; i < MAX_DEPTH; ++i) {
		DownsampleGraph(mAdj[i], mV[i], mN[i], mA[i], mV[i + 1], mN[i + 1], mA[i + 1],
			mToUpper[i], mToLower[i], mAdj[i + 1]);
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
	mScale = scale;
	for (size_t i = 0; i<mV.size(); ++i) {
		mQ[i].resize(mN[i].rows(), mN[i].cols());
		mO[i].resize(mN[i].rows(), mN[i].cols());
		for (int j = 0; j < mN[i].cols(); ++j) {
			Vector3f s, t;
			coordinate_system(mN[i].col(j), s, t);
			float angle = ((float)rand()) / RAND_MAX * 2 * M_PI;
			float x = ((float)rand()) / RAND_MAX * 2 - 1.f;
			float y = ((float)rand()) / RAND_MAX * 2 - 1.f;
			mQ[i].col(j) = s * std::cos(angle) + t * std::sin(angle);
			mO[i].col(j) = mV[i].col(j) + (s * x + t * y) * scale;
		}
	}
}

void Hierarchy::DownsampleGraph(const AdjacentMatrix adj, const MatrixXf &V,
	const MatrixXf &N, const VectorXf &A,
	MatrixXf &V_p, MatrixXf &N_p, VectorXf &A_p,
	MatrixXi &to_upper, VectorXi &to_lower,
	AdjacentMatrix& adj_p) {
	struct Entry {
		int i, j;
		float order;
		inline Entry() { i = j = -1; };
		inline Entry(int i, int j, float order) : i(i), j(j), order(order) { }
		inline bool operator<(const Entry &e) const { return order > e.order; }
		inline bool operator==(const Entry &e) const { return order == e.order; }
	};
	int nLinks = 0;
	for (auto& adj_i : adj)
		nLinks += adj_i.size();
	std::vector<Entry> entries(nLinks);
	int offset = 0;
	for (int i = 0; i < V.cols(); ++i) {
		for (int j = 0; j < adj[i].size(); ++j) {
			int k = adj[i][j].id;
			float dp = N.col(i).dot(N.col(k));
			float ratio = A[i] > A[k] ? (A[i] / A[k]) : (A[k] / A[i]);
			entries[offset + j] = Entry(i, k, dp * ratio);
		}
		offset += adj[i].size();
	}

	std::stable_sort(entries.begin(), entries.end());
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

	for (int i = 0; i < nCollapsed; ++i) {
		const Entry &e = entries[i];
		const float area1 = A[e.i], area2 = A[e.j], surfaceArea = area1 + area2;
		if (surfaceArea > RCPOVERFLOW)
			V_p.col(i) = (V.col(e.i) * area1 + V.col(e.j) * area2) / surfaceArea;
		else
			V_p.col(i) = (V.col(e.i) + V.col(e.j)) * 0.5f;
		Vector3f normal = N.col(e.i) * area1 + N.col(e.j) * area2;
		float norm = normal.norm();
		N_p.col(i) = norm > RCPOVERFLOW ? Vector3f(normal / norm)
			: Vector3f::UnitX();
		A_p[i] = surfaceArea;
		to_upper.col(i) << e.i, e.j;
		to_lower[e.i] = i; to_lower[e.j] = i;
	}

	offset = nCollapsed;

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
	for (int i = 0; i < V_p.cols(); ++i) {
		std::vector<Link> scratch;
		for (int j = 0; j<2; ++j) {
			int upper = to_upper(j, i);
			if (upper == -1)
				continue;
			for (auto& link : adj[upper])
				scratch.push_back(Link(to_lower[link.id], link.weight));
		}
		std::sort(scratch.begin(), scratch.end());
		int id = -1;
		for (auto& link : scratch) {
			if (link.id != i) {
				if (id != link.id) {
					adj_p[i].push_back(link);
					id = link.id;
				}
				else {
					adj_p[i].back().weight += link.weight;
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
}
