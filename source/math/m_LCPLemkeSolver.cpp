#include "m_LCPLemkeSolver.h"

#include "m_Const.h"

#include "helpers/debug.h"

using namespace math;
/*__int64 GetCPUClock()
{
    __int64 res;
    __asm
	{
        rdtsc
        mov dword ptr res, eax
        mov dword ptr res+4, edx
    }
    return res;
}*/
CLCPLemkeSolver::CLCPLemkeSolver()
{
	mSize = 0;
	mpInvbasmat = NULL;
	mpBasis = NULL;
	mpResult = NULL;
	mpTransformedC = NULL;
	mpTransformedQ = NULL;
	mpIndices0 = NULL;
	mpIndices1 = NULL;

	return;
}
CLCPLemkeSolver::~CLCPLemkeSolver()
{
	Empty();
	return;
}

bool CLCPLemkeSolver::SetSize(const unsigned long size)
{
	if(size <= mSize)
			{
		mSize = size;
		}
	else
			{
		this->Empty();
		mSize = size;
		mpInvbasmat = new FLOATTYPE [mSize*mSize];
		mpBasis = new SLemkeVar[mSize];
		mpResult = new FLOATTYPE* [mSize]; 
		mpTransformedC = new FLOATTYPE[mSize];
		mpTransformedQ = new FLOATTYPE[mSize];
		mpIndices0 = new unsigned long[mSize];
		mpIndices1 = new unsigned long[mSize];
		}

	return true;
}
void CLCPLemkeSolver::Empty(void)
{
	if(mpInvbasmat != NULL) 
		delete[] mpInvbasmat;
	if(mpBasis != NULL)
		delete[] mpBasis;
	if(mpResult != NULL)
		delete[] mpResult;
	if(mpTransformedC != NULL)
		delete[] mpTransformedC;
	if(mpTransformedQ != NULL)
		delete[] mpTransformedQ;
	if(mpIndices0 != NULL)
		delete[] mpIndices0;
	if(mpIndices1 != NULL)
		delete[] mpIndices1;

	mSize = 0;
	return;
}
unsigned long CLCPLemkeSolver::GetSize(void)
{
	return mSize;
}

bool CLCPLemkeSolver::FastInverseMatrix(unsigned long new_var_pos)
{
	unsigned long j;

	mpTransformedC[new_var_pos] -= 1.0;
	FLOATTYPE denom = mpTransformedC[new_var_pos] + 1.0f;
	if(denom == 0)
	{
		return false;
	}
	denom = 1.0f / denom;
	for(unsigned long i = 0; i < new_var_pos; ++i)
	{
		for( j = 0; j < mSize; ++j)
		{
			mpInvbasmat[i * mSize + j] -= (denom * mpTransformedC[i] * mpInvbasmat[new_var_pos * mSize + j]);
		}
	}
	for(unsigned long i = new_var_pos + 1; i < mSize; ++i)
	{

		for( j = 0; j < mSize; ++j)
		{
			mpInvbasmat[i * mSize + j] -= (denom * mpTransformedC[i] * mpInvbasmat[new_var_pos * mSize + j]);
		}
	}
	for( j = 0; j < mSize; ++j)
	{
		mpInvbasmat[new_var_pos * mSize + j] -= (denom * mpTransformedC[new_var_pos] * mpInvbasmat[new_var_pos * mSize + j]);
	}
	return true;
}

bool CLCPLemkeSolver::FindDeparting(CMatrixN &LCPMat, CVectorN &LCPVec, SLemkeVar entering_varible, unsigned long &departing_varible)
{
	unsigned long i, j;
	ZeroMemory(mpTransformedC, mSize * sizeof(FLOATTYPE));
	if(entering_varible.type == LEMKE_VAR_Y)
	{
		for(i = 0; i < mSize; ++i)
		{
			for(j = 0;j < mSize; ++j)
			{
				mpTransformedC[i] -= mpInvbasmat[i * mSize +j] * LCPMat.Get(j, entering_varible.index);
			}
		}
	}
	else if(entering_varible.type == LEMKE_VAR_W)
	{
		for(j = 0; j < mSize; ++j)
		{
			mpTransformedC[j] = mpInvbasmat[j * mSize + entering_varible.index];
		}
	}
	else
	{
		assert(false && "Wrong varible type in Find Departing");
		return false;
	}
	ZeroMemory(mpTransformedQ, mSize * sizeof(FLOATTYPE));
	for(i = 0;i < mSize; ++i)
	{
		for(j = 0;j < mSize; ++j)
		{
			mpTransformedQ[i] += mpInvbasmat[i * mSize +j] * LCPVec.Get(j);
		}
	}	
	unsigned long indsize = 0, indsize1;
	unsigned long *srcindsize,*dstindsize ;
	FLOATTYPE minimum, curr;
	unsigned long **pIndSrc, **pIndDest;
	pIndDest = &mpIndices0;
	bool FoundAny = false;
	for( unsigned long i = 0; i < mSize; i++)
	{
		if(mpTransformedC[i] > 0)
		{
			curr = mpTransformedQ[i] / mpTransformedC[i];
			if((!FoundAny) || (FoundAny && ( M_EQUAL_AT(curr , minimum, 0.01*M_EPSILON) ) ) )
			{
				FoundAny = true;
				(*pIndDest)[indsize] = i;
				++indsize;
				minimum = curr;
			}
			else if(FoundAny && ( M_LESS_AT(curr , minimum, 0.01 * M_EPSILON) ))
			{
				minimum = curr;
				(*pIndDest)[0] = i;
				indsize = 1;
			}
		}
	}
	if(indsize == 0)
	{
		return false;
	}
	else if(indsize == 1)
	{
		departing_varible = (*pIndDest)[0];
		return true;
	}
	for(unsigned long i = 0 ; i < mSize; i++)
	{
		//for each column of InvMat
		if(i % 2 == 0)
		{
			pIndSrc = &mpIndices0;
			pIndDest = &mpIndices1;
			srcindsize = &indsize;
			dstindsize = &indsize1;
		
		}
		else
		{
			pIndSrc = &mpIndices1;
			pIndDest = &mpIndices0;
			srcindsize = &indsize1;
			dstindsize = &indsize;
		
		}
		(*dstindsize) = 0;
		FoundAny = false;
		minimum = mpInvbasmat[(*pIndSrc)[0] * mSize + i] / mpTransformedC[(*pIndSrc)[0]];
		for(unsigned long j = 0; j < (*srcindsize); j++)
		{
			if(mpBasis[(*pIndSrc)[j]].type == LEMKE_VAR_Z0)
			{
				departing_varible = (*pIndSrc)[j];
				return true;
			}
			curr = mpInvbasmat[(*pIndSrc)[j] * mSize +i] / mpTransformedC[(*pIndSrc)[j]];
			//Warning! "Smart" changes made late at night
			if(FoundAny && M_EQUAL_AT(curr, minimum, 0.01 * M_EPSILON))
			{
				(*pIndDest)[(*dstindsize)] = (*pIndSrc)[j];
				++(*dstindsize);
			}
			else if(!FoundAny || (FoundAny && M_LESS_AT(curr, minimum, 0.01 * M_EPSILON)))
			{
				FoundAny = true;
				minimum = curr;
				(*pIndDest)[0] = (*pIndSrc)[j];
				(*dstindsize) = 1;
			}
		}
		if((*dstindsize) == 1) break;
	}
	if((*dstindsize) != 1)
	{
	//	assert(false && "Couldnt remove ambigiousity in FindDeparting");
	}
	departing_varible = (*pIndDest)[0];
	return true;
}

CLCPLemkeSolver::ELCPErr CLCPLemkeSolver::Solve(CMatrixN &LCPMat, CVectorN &LCPVec, CVectorN &y, CVectorN &w)
{
//	__int64 Time0 = GetCPUClock();;
	unsigned long start_index = 0;
	FLOATTYPE minimum = LCPVec.Get(0);
	mpBasis[0].index = 0;
	mpBasis[0].type = LEMKE_VAR_W;

	for(unsigned long i = 1; i < mSize; i++)
	{
		if(LCPVec.Get(i) < minimum)
		{
			minimum = LCPVec.Get(i);
			start_index = i;
		}
		mpBasis[i].index = i;
		mpBasis[i].type = LEMKE_VAR_W;
	}
	if(minimum >= 0)
	{
		y.MakeZero();
		w = LCPVec;
		return LCPERR_OK;
	}

	mpBasis[start_index].type = LEMKE_VAR_Z0;

	ZeroMemory(mpInvbasmat, mSize * mSize * sizeof(FLOATTYPE));
	for(unsigned long j = 0; j < mSize; j++)
	{
		mpInvbasmat[j * mSize + j] = 1.0;
		mpInvbasmat[j * mSize + start_index] = -1.0;
	}

	SLemkeVar entering_varible, temp_varible;
	unsigned long departing_varible;
	entering_varible.index = start_index;
	entering_varible.type = LEMKE_VAR_Y;
	unsigned long cycle_count = 0;
	while(true)
	{
		cycle_count++ ;
		if(cycle_count > 1000)
		{
/*DEBUG*/	assert("Overcycling in Lemke" && false);
			return LCPERR_OVERCYCLING;
		}
	//	__int64 Time1 = GetCPUClock();
		if(!FindDeparting(LCPMat, LCPVec, entering_varible, departing_varible))
		{
/*DEBUG*/	assert("Lemeke failed!" && false);
			printf("\n\n");
			for(unsigned long i = 0; i < mSize; i++)
			{
					printf("%f ", mpTransformedC[i]);
			}
			printf("\n\n");
			for(unsigned long i = 0; i < mSize; i++)
			{
					printf("%f ", mpTransformedQ[i]);
			}
			printf("\n\n");
			for(unsigned long i = 0; i < mSize; i++)
			{
					printf("%f ", mpTransformedQ[i] / mpTransformedC[i]);
			}
			return LCPERR_CANTSOLVE;
		}

		if(!FastInverseMatrix(departing_varible))
		{
/*DEBUG*/	assert(false && "Basis unfeasable!");
			return LCPERR_CANTINVERT;
		}

		bool bEnd = mpBasis[departing_varible].type == LEMKE_VAR_Z0;
		temp_varible = entering_varible;
		entering_varible.index = mpBasis[departing_varible].index;
		entering_varible.type = !(mpBasis[departing_varible].type);

		mpBasis[departing_varible] = temp_varible;
		if(bEnd)
		{
			break;
		}
	}
	w.MakeZero();
	y.MakeZero();
	
	FLOATTYPE temp;
	for(unsigned long i = 0; i < mSize; ++i)
	{
		temp = 0.0;
		for(unsigned long j = 0; j < mSize; ++j)
		{
			temp += mpInvbasmat[i * mSize +j] * LCPVec.Get(j);
		}
		if(mpBasis[i].type == LEMKE_VAR_Y)
		{
			y[mpBasis[i].index] = temp;
		}
		else
		{
			w[mpBasis[i].index] = temp;
		}
	}

	return LCPERR_OK;
}

