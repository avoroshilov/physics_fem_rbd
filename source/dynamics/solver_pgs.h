#ifndef _DYNAMICS_SOLVER_PGS_H_
#define _DYNAMICS_SOLVER_PGS_H_

#include "solver_base.h"

class SolverPGS : public SolverBase
{
	friend class FEMJoint;		// list all joints here; in future make some kind of interface for SOlver tompass to joints
	friend class BallSocket;	// list all joints here; in future make some kind of interface for SOlver tompass to joints

public:

	SolverPGS(unsigned int MaxIterations, float Precision) : SolverBase(MaxIterations, Precision)
	{
	}

	void Solve(float dt);

	dynarray<float, MAX_NUM_JACOBI> m_invDiag;							// D^-1 from L + D + U
};

#endif