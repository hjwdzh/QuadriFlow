#ifndef OPTIMIZER_H_
#define OPTIMIZER_H_
#include "hierarchy.h"

class Optimizer
{
public:
	Optimizer();
	static void optimize_orientations(Hierarchy &mRes);
	static void optimize_scale(Hierarchy &mRes);
	static void optimize_positions(Hierarchy &mRes);
};
#endif