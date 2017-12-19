#ifndef PARAMETRIZER_H_
#define PARAMETRIZER_H_

#include <Eigen/Core>
#include <Eigen/Dense>
#include <map>
#include "AdjacentMatrix.h"
#include "hierarchy.h"

using namespace Eigen;

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

	// Extract Mesh
	void ExtractMesh();

	std::map<int, int> vertex_singularities;
	std::map<int, int> singularities;
	// input mesh
	MatrixXf V;
	MatrixXf N;
	MatrixXi F;

	// data structures
	VectorXi V2E;
	VectorXi E2E;
	VectorXi boundary;
	VectorXi nonManifold;
	AdjacentMatrix adj;
	Hierarchy hierarchy;
	
	// Mesh Status;
	float surface_area;
	float scale;
	float average_edge_length;
	float max_edge_length;
	VectorXf A;

	// target mesh
	int num_vertices;
	int num_faces;

	// Extracted Mesh
	std::vector<std::vector<TaggedLink>> adj_extracted;

	MatrixXi mF_extracted;
	MatrixXf mV_extracted;
	MatrixXf mN_extracted, mNf_extracted;

	std::vector<VectorXi> qF;
	std::vector<std::vector<std::pair<int, int> > > qVF;
};
#endif