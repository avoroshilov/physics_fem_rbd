#ifndef _DYNAMICS_SOLVER_CG_H_
#define _DYNAMICS_SOLVER_CG_H_

#include "solver_base.h"

class SolverCG : public SolverBase
{
	friend class FEMJoint;		// list all joints here; in future make some kind of interface for SOlver tompass to joints
	friend class BallSocket;	// list all joints here; in future make some kind of interface for SOlver tompass to joints

public:

	SolverCG(unsigned int MaxIterations, float Precision) : SolverBase(MaxIterations, Precision)
	{
	}

	void Solve(float dt);

	dynarray<float, MAX_NUM_JACOBI> m_ad_vec;								// vector A * d_i
	dynarray<float, MAX_NUM_JACOBI> m_r;									// residual
	dynarray<float, MAX_NUM_JACOBI> m_d;									// direction
};

#endif