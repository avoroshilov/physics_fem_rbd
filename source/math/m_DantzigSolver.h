#ifndef _M_DANTZIG_SOLVER_H_
#define _M_DANTZIG_SOLVER_H_

#include "all.h"

#include "m_MatrixN.h"
#include "m_SetOfIndices.h"
#include "m_VectorN.h"

namespace math
{
	class CDantzigSolver
	{
	public:
		CMatrixN							mAMatrix;
		CVectorN							mBVector;
		CVectorN							mResult;

		explicit CDantzigSolver(unsigned long size);
		~CDantzigSolver(void);

		bool								Solve(void);

	private:
		bool								DriveToZero(unsigned long i);
		bool								ComputeStepDirection(unsigned long i);
		bool								ComputeMaxStep(unsigned long i, FLOATTYPE &s, unsigned long &j);

	private:
		unsigned long						mSize;

		CVectorN							mAVector;
		CSetOfIndices						mSetC;
		CSetOfIndices						mSetNC;

		CVectorN							mDeltaResult;
		CVectorN							mDeltaAVector;
	};
}

#endif
