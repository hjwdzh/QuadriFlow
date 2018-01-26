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
	void Initialize(double scale);
	void DownsampleGraph(const AdjacentMatrix adj, const MatrixXd &V,
		const MatrixXd &N, const VectorXd &A,
		MatrixXd &V_p, MatrixXd &N_p, VectorXd &A_p,
		MatrixXi &to_upper, VectorXi &to_lower,
		AdjacentMatrix& adj_p);
	void generate_graph_coloring_deterministic(const AdjacentMatrix &adj, int size,
		std::vector<std::vector<int> > &phases);

	enum {
		MAX_DEPTH = 25
	};

	void SaveToFile(FILE* fp);
	void LoadFromFile(FILE* fp);

	double mScale;

	MatrixXi mF;
	VectorXi mE2E;
	std::vector<AdjacentMatrix> mAdj;
	std::vector<MatrixXd> mV;
	std::vector<MatrixXd> mN;
	std::vector<VectorXd> mA;
	std::vector<VectorXi> mToLower;
	std::vector<MatrixXi> mToUpper;
	std::vector<std::vector<std::vector<int> > > mPhases;
	// parameters
	std::vector<MatrixXd> mQ;
	std::vector<MatrixXd> mO;
	std::vector<MatrixXd> mS;
	std::vector<MatrixXd> mK;
};
#endif