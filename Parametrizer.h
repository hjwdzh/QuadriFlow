#ifndef PARAMETRIZER_H_
#define PARAMETRIZER_H_
#include <list>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <map>
#include <set>
#include "AdjacentMatrix.h"
#include "hierarchy.h"
#include "serialize.h"
#include "DisajointTree.h"
using namespace Eigen;

typedef std::pair<unsigned int, unsigned int> Edge;
typedef std::map<int, std::pair<int, int> > SingDictionary;

struct DEdge
{
	DEdge()
		: x(0), y(0)
	{}
	DEdge(int _x, int _y) {
		if (_x > _y)
			x = _y, y = _x;
		else
			x = _x, y = _y;
	}
	bool operator<(const DEdge& e) const {
		return x < e.x || x == e.x && y < e.y;
	}
	bool operator==(const DEdge& e) const {
		return x == e.x && y == e.y;
	}
	bool operator!=(const DEdge& e) const {
		return x != e.x || y != e.y;
	}
	int x, y;
};
struct ExpandInfo
{
	ExpandInfo()
	{}
	int current_v;
	int singularity;
	int step;
	int edge_id;
	int prev;
};

class Parametrizer
{
public:
	void Load(const char* filename);
	void Initialize();
	
	// member function
	void ComputeMeshStatus();
	void ComputeSmoothNormal();
	void ComputeVertexArea();
	void ComputeOrientationSingularities();
	void ComputePositionSingularities(int with_scale = 0);
	void EstimateScale();

	// index map
	void ComputeIndexMap(int with_scale = 0);
	void BuildEdgeInfo();
	void ComputeMaxFlow();
	void BuildIntegerConstraints();
	void ComputePosition(int with_scale = 0);
	void FixFlipAdvance();

	// sanity check
	void SanityCheckDiff(int sing);

	// File IO
	void SaveToFile(FILE* fp);
	void LoadFromFile(FILE* fp);
	void ExtractMesh(const char* obj_name);

	std::map<int, int> singularities;
	std::map<int, Vector2i> pos_sing;
	MatrixXi pos_rank;
	MatrixXi pos_index;
	// input mesh
	MatrixXd V;
	MatrixXd N;
	MatrixXd Nf;
	MatrixXd FS;
	MatrixXd FQ;
	MatrixXi F;

	std::vector<MatrixXd> triangle_space;

	// data structures
	VectorXi V2E;
	VectorXi E2E;
	VectorXi boundary;
	VectorXi nonManifold;
	AdjacentMatrix adj;
	Hierarchy hierarchy;
	
	// Mesh Status;
	double surface_area;
	double scale;
	double average_edge_length;
	double max_edge_length;
	VectorXd A;

	// target mesh
	int num_vertices;
	int num_faces;

	//just for test
	DisajointTree disajoint_tree;
	DisajointOrientTree disajoint_orient_tree;

	int compact_num_v;
	std::vector<Vector3d> O_compact;
	std::vector<Vector4i> F_compact;
	std::vector<int> bad_vertices;
	std::vector<double> counter;
	std::vector<Vector2i> edge_diff;
	std::map<DEdge, int> edge_ids;
	std::vector<DEdge> edge_values;

	std::vector<Vector3i> constraints_index;
	std::vector<Vector3i> constraints_sign;
	std::vector<std::pair<Vector2i, int> > variables;
	std::vector<Vector4i> edge_to_constraints;

	std::set<DEdge> cuts;
	std::vector<Vector3i> flipped;

	// fixed_vertices
	std::vector<int> fixed;
	std::set<DEdge> fixed_cuts;
};
#endif