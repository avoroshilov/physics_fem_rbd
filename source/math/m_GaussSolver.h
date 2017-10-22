#ifndef _GAUSS_SOLVER_H_
#define _GAUSS_SOLVER_H_

#include "all.h"

#include "m_MatrixN.h"

namespace math
{
	class CGaussSolver
	{
	public:
		CMatrixN							mMatrix;
		CVectorN							mVector;
		CVectorN							mResult;

		CGaussSolver(void){};
		explicit CGaussSolver(unsigned long n) : mMatrix(n), mVector(n), mResult(n){};
		~CGaussSolver(void){};

		bool								SetSize(unsigned long size);

		bool								Solve(void);
		bool								Solve0(void);
		bool								Solve1(void);
	};
}

#endif
