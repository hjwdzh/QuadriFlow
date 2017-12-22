#ifndef PARAMETRIZER_H_
#define PARAMETRIZER_H_

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
	void ExtractMesh();
	void LoopFace(int mode);

	void SaveToFile(FILE* fp);

	void LoadFromFile(FILE* fp);

	std::map<int, int> vertex_singularities;
	std::map<int, int> singularities;
	
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
};
#endif