#ifndef _DYNAMICS_SOLVER_LCPCG_H_
#define _DYNAMICS_SOLVER_LCPCG_H_

#include "solver_base.h"

class SolverLCPCG : public SolverBase
{
	friend class FEMJoint;		// list all joints here; in future make some kind of interface for SOlver tompass to joints
	friend class BallSocket;	// list all joints here; in future make some kind of interface for SOlver tompass to joints

public:

	SolverLCPCG(unsigned int MaxIterations, float Precision) : SolverBase(MaxIterations, Precision)
	{
	}

	void Solve(float dt);

	//dynarray<float, MAX_NUM_JACOBI> m_invDiag;							// inverse diagonal

	dynarray<float, MAX_NUM_JACOBI> m_ad_vec;								// vector A * d
	dynarray<float, MAX_NUM_JACOBI> m_ap_vec;								// vector A * p

	dynarray<float, MAX_NUM_JACOBI> m_p;									// stepping vector
	dynarray<float, MAX_NUM_JACOBI> m_d;									// stepping vector, also intermediate in Power Iteration
	dynarray<float, MAX_NUM_JACOBI> m_g;									// gradient vector
};

#endif