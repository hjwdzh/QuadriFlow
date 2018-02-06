#ifndef OPTIMIZER_H_
#define OPTIMIZER_H_
#include "config.hpp"
#include "hierarchy.hpp"

class Optimizer
{
public:
	Optimizer();
	static void optimize_orientations(Hierarchy &mRes);
	static void optimize_scale(Hierarchy &mRes);
	static void optimize_positions(Hierarchy &mRes, int with_scale = 0);
    static void optimize_integer_constraints(Hierarchy &mRes, std::map<int, int>& singularities);
    
#ifdef WITH_CUDA
	static void optimize_orientations_cuda(Hierarchy &mRes);
	static void optimize_positions_cuda(Hierarchy &mRes);
#endif
};

#ifdef WITH_CUDA
extern void UpdateOrientation(int* phase, int num_phases, glm::dvec3* N, glm::dvec3* Q, Link* adj, int* adjOffset, int num_adj);
extern void PropagateOrientationUpper(glm::dvec3* srcField, int num_orientation, glm::ivec2* toUpper, glm::dvec3* N, glm::dvec3* destField);
extern void PropagateOrientationLower(glm::ivec2* toUpper, glm::dvec3* Q, glm::dvec3* N, glm::dvec3* Q_next, glm::dvec3* N_next, int num_toUpper);

extern void UpdatePosition(int* phase, int num_phases, glm::dvec3* N, glm::dvec3* Q, Link* adj, int* adjOffset, int num_adj, glm::dvec3* V, glm::dvec3* O, double scale);
extern void PropagatePositionUpper(glm::dvec3* srcField, int num_position, glm::ivec2* toUpper, glm::dvec3* N, glm::dvec3* V, glm::dvec3* destField);

extern void JacobiSolve(std::vector<double>& D, std::vector<double>& R, std::vector<int>& R_ind, std::vector<int>& R_offset, std::vector<double>& x, std::vector<double>& b);
#endif

#endif
