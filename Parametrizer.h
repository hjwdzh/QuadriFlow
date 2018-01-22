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
	// Extract Mesh
	void ExtractMesh(int with_scale = 0);
	void LoopFace(int mode);

	void SaveToFile(FILE* fp);

	void LoadFromFile(FILE* fp);

	std::map<int, int> vertex_singularities;
	std::map<int, int> singularities;
	std::map<int, Vector2i> pos_sing;
	MatrixXi pos_rank;
	MatrixXi pos_index;
	std::vector<DEdge> valid_edges;
	std::vector<DEdge> mE_extracted, mE_extracted2;
	std::vector<Vector4i> mF_extracted2;
	std::vector<Vector3d> mO_extracted;
	std::set<int> color_tests[4];
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

	// Extracted Mesh
	std::vector<std::vector<TaggedLink>> adj_extracted;

	MatrixXi mF_extracted;
	MatrixXd mV_extracted;
	MatrixXd mN_extracted, mNf_extracted;

	std::vector<VectorXi> qF;
	std::vector<std::vector<std::pair<int, int> > > qVF;

	// singularity graph
	std::map<Edge, int> edge_idmap;
	std::vector<Edge> qE;
	std::vector<std::vector<int> > qVE;
	std::vector<std::vector<int> > qEV;
	std::vector<int> qEE, qRE;
	std::vector<SingDictionary> sin_graph;
	std::vector<std::set<int> > triangle_edge_pair;
	std::vector<ExpandInfo > q;
	int front_index;

	// singularity
	std::vector<std::vector<int> > edge_strips;
	std::vector<std::set<std::pair<int, int> > > singularity_entry;

	// index map
	void ComputeIndexMap(int with_scale = 0);
	std::vector<int> vertex_rank;

	//just for test
	std::vector<int> colors;
	std::vector<Vector3d> shrink_pts;
	std::set<int> singular_patches;
	std::set<DEdge> singular_e;

	int compact_num_v;
	std::map<Edge, std::pair<int, Vector2i> > Es;
	std::vector<std::pair<int, int> > shrink_parents;
	std::vector<int> shrink_compact_indices;
	std::vector<int> shrink_ranks;
	std::vector<Vector3d> O_compact;
	std::vector<double> counter;
	std::vector<Vector2i> edge_diff;
	std::map<DEdge, int> edge_ids;
	std::set<int> singular_patches_buf;
	std::set<DEdge> singular_e_buf;
	std::vector<DEdge> edge_values;
	std::set<DEdge> singular_edges;
	std::vector<std::list<int> > edge_neighbors;

	std::vector<Vector3i> constraints_index;
	std::vector<Vector3i> constraints_sign;
	std::vector<std::pair<Vector2i, int> > variables;
	std::vector<std::pair<int, int> > parentss;
	std::vector<Vector4i> edge_to_constraints;

	std::vector<Vector2i> param;
	Vector2i min_param, max_param;
	void ComputeMaxFlow();
	void MergeVertices(int v);
	void BuildIntegerConstraints();
	void UpdateMesh();
	void FixFlip();
	std::set<DEdge> cuts;

	std::vector<int> vdist;

};
#endif