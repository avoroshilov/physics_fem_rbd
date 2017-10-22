#ifndef _M_LCPLEMKESOLVER_H_
#define _M_LCPLEMKESOLVER_H_

#include "all.h"

#include "m_MatrixN.h"
#include "m_VectorN.h"

namespace math
{
	const unsigned char LEMKE_VAR_Y(0x0), LEMKE_VAR_W(0x1), LEMKE_VAR_Z0(0x2);

	class CLCPLemkeSolver
	{
	public:
		
		CLCPLemkeSolver();
		~CLCPLemkeSolver();
		enum ELCPErr {LCPERR_OK = 0, LCPERR_OVERCYCLING = 1, LCPERR_CANTSOLVE = 2, LCPERR_CANTINVERT = 3, LCPERR_CANTSOLVELES = 4}; 
		
		ELCPErr Solve(CMatrixN &LCPMat, CVectorN &LCPVec, CVectorN &y, CVectorN &w);

		bool SetSize(const unsigned long size);
		void Empty(void);
		unsigned long GetSize(void); 
	private:
		struct SLemkeVar
		{
			unsigned long index;
			unsigned char type:2;
		};
		unsigned long mSize;
		FLOATTYPE *mpInvbasmat;
		SLemkeVar *mpBasis;
		FLOATTYPE **mpResult;
		FLOATTYPE *mpTransformedC;
		FLOATTYPE *mpTransformedQ;
		unsigned long *mpIndices0;
		unsigned long *mpIndices1;

		bool FastInverseMatrix(unsigned long new_var_pos);
		bool FindDeparting(CMatrixN &LCPMat, CVectorN &LCPVec, SLemkeVar entering_varible, unsigned long &departing_varible);
	};
}

#endif
