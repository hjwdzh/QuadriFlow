#ifndef HIERARCHY_H_
#define HIERARCHY_H_

#include <vector>

#include "AdjacentMatrix.h"
#include "serialize.h"
#define RCPOVERFLOW   2.93873587705571876e-39f

class Hierarchy
{
public:
	Hierarchy();
	void Initialize(float scale);
	void DownsampleGraph(const AdjacentMatrix adj, const MatrixXf &V,
		const MatrixXf &N, const VectorXf &A,
		MatrixXf &V_p, MatrixXf &N_p, VectorXf &A_p,
		MatrixXi &to_upper, VectorXi &to_lower,
		AdjacentMatrix& adj_p);

	enum {
		MAX_DEPTH = 25
	};

	void SaveToFile(FILE* fp);
	void LoadFromFile(FILE* fp);

	float mScale;

	MatrixXi mF;
	VectorXi mE2E;
	std::vector<AdjacentMatrix> mAdj;
	std::vector<MatrixXf> mV;
	std::vector<MatrixXf> mN;
	std::vector<VectorXf> mA;
	std::vector<VectorXi> mToLower;
	std::vector<MatrixXi> mToUpper;

	// parameters
	std::vector<MatrixXf> mQ;
	std::vector<MatrixXf> mO;
	std::vector<MatrixXf> mS;
};
#endif