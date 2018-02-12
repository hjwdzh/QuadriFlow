#ifndef PARAMETRIZER_H_
#define PARAMETRIZER_H_
#include <atomic>
#include <condition_variable>
#ifdef WITH_TBB
#include <tbb/tbb.h>
#endif

#include <list>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <map>
#include <set>
#include "adjacent-matrix.hpp"
#include "field-math.hpp"
#include "hierarchy.hpp"
#include "post-solver.hpp"
#include "serialize.hpp"
#include "disajoint-tree.hpp"
using namespace Eigen;

typedef std::pair<unsigned int, unsigned int> Edge;
typedef std::map<int, std::pair<int, int> > SingDictionary;

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
    // Mesh Initialization
	void Load(const char* filename);
    void ComputeMeshStatus();
    void ComputeSmoothNormal();
    void ComputeVertexArea();
	void Initialize(int faces, int with_scale = 0);
	
	// Singularity and Mesh property
	void ComputeOrientationSingularities();
	void ComputePositionSingularities(int with_scale = 0);
#ifdef WITH_SCALE
	void EstimateScale();
#endif
    
	// Integer Grid Map Pipeline
	void ComputeIndexMap(int with_scale = 0);
	void BuildEdgeInfo();
	void ComputeMaxFlow();
	void BuildIntegerConstraints();

    // Fix Flip
    void FixFlipHierarchy();
    void FixHoles();
    double QuadEnergy(std::vector<int>& loop_vertices, std::vector<Vector4i>& res_quads, int level);

    // Quadmesh and IO
    void ExtractQuadMesh();
	void OutputMesh(const char* obj_name);

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

	int compact_num_v;
	std::vector<Vector3d> O_compact;
    std::vector<Vector3d> Q_compact;
    std::vector<Vector3d> N_compact;
	std::vector<Vector4i> F_compact;
    VectorXi V2E_compact;
    std::vector<int> E2E_compact;
    VectorXi boundary_compact;
    VectorXi nonManifold_compact;
    
	std::vector<int> bad_vertices;
	std::vector<double> counter;
	std::vector<Vector2i> edge_diff;
	std::vector<DEdge> edge_values;
	std::vector<Vector3i> face_edgeIds;
	std::vector<Vector3i> face_edgeOrients;

	std::vector<std::pair<Vector2i, int> > variables;

	std::set<DEdge> cuts;
	std::vector<Vector3i> flipped;

	// fixed_vertices
	std::vector<int> fixed;
	std::set<DEdge> fixed_cuts;

	std::set<int> edge_around_singularities;

};

extern void generate_adjacency_matrix_uniform(const MatrixXi& F, const VectorXi& V2E,
                                              const VectorXi& E2E, const VectorXi& nonManifold,
                                              AdjacentMatrix& adj);


#endif
