#include "m_DantzigSolver.h"

#include "m_Const.h"
#include "m_GaussSolver.h"

#include <stdio.h>

#include "helpers/log.h"

#pragma message(__FILE__" : warning: Rewrite code using M_EPSILON!")

using namespace math;

CDantzigSolver::CDantzigSolver(unsigned long size) : mSetC(size), mSetNC(size)
{
	mSize = size;

	mAMatrix.SetSize(size);
	mBVector.SetSize(size);
	mResult.SetSize(size);

	mAVector.SetSize(size);
	mDeltaResult.SetSize(size);
	mDeltaAVector.SetSize(size);
}

CDantzigSolver::~CDantzigSolver(void)
{
}

bool CDantzigSolver::Solve(void)
{
	gLog.Print("--------------------------------------------------------------------------------");

	mAMatrix.Log();

	gLog.Print("");

	mBVector.Log();

	gLog.Print("--------------------------------------------------------------------------------");

	mResult.MakeZero();
	mAVector = mBVector;

	//printf("%.5f %.5f  %.5f %.5f", mAVector[0], mAVector[1], mBVector[0], mBVector[1]);

	mSetC.Empty();
	mSetNC.Empty();

	for(unsigned long i = 0; i < mSize; i++)
	{
		gLog.Print("Stage %ld", i);
		//printf("VecA size: %ld, VecB: %ld", mAVector.GetSize(), mBVector.GetSize());
		//printf("mAVector[i] = %.5f", mAVector[i]);

		gLog.Print("Vector a:");
		mAVector.Log();

		if(mAVector[i] < 0)
		{
			gLog.Print("i = %ld pass.", i);

			if(!DriveToZero(i))
			{
				gLog.Print("DriveToZero() failed.");
				gLog.Print("FFailed to solve");
				mAMatrix.Log();
				gLog.Print("");
				mBVector.Log();
				mResult.Log();

				return false;
			}
		}
	}

	gLog.Print("Solved OK");
	mResult.Log();

	return true;
}

bool CDantzigSolver::DriveToZero(unsigned long i)
{
	unsigned long counter = 0;

	while(true)
	{
		if(counter > 20)
		{
			gLog.Print("Overcycling.");

			mAMatrix.Log();
			gLog.Print("");
			mAVector.Log();
			mResult.Log();

			gLog.Print("mSetC:");

			for(unsigned long n = 0; n < mSetC.GetSize(); n++)
				gLog.Print("%ld", mSetC[n]);

			gLog.Print("");
			gLog.Print("mSetNC:");

			for(unsigned long n = 0; n < mSetNC.GetSize(); n++)
				gLog.Print("%ld", mSetNC[n]);

			return false;
		}

		gLog.Print("----------Step: %ld------------------------------", counter);
		counter++;

		gLog.Print("Matrix_a:");
		mAMatrix.Log();
		gLog.Print("");

		gLog.Print("Vector_a:");
		mAVector.Log();
		gLog.Print("");

		gLog.Print("Vector_b:");
		mBVector.Log();
		gLog.Print("");

		gLog.Print("Result:");
		mResult.Log();
		gLog.Print("");


		gLog.Print("mSetC:");

			for(unsigned long n = 0; n < mSetC.GetSize(); n++)
				gLog.Print("%ld", mSetC[n]);

			gLog.Print("");
			gLog.Print("mSetNC:");

			for(unsigned long n = 0; n < mSetNC.GetSize(); n++)
				gLog.Print("%ld", mSetNC[n]);

		if(!ComputeStepDirection(i))
		{
			gLog.Print("Compute direction failed.");
	
			return false;
		}

		gLog.Print("Delta r vector:");
		mDeltaResult.Log();


		mDeltaAVector = mAMatrix * mDeltaResult;

		gLog.Print("Delta a vector:");
		mDeltaAVector.Log();

		//printf("DeltaA = %.5f", mDeltaAVector[0]);
		//printf("AA = %.5f", mAMatrix(0, 0));

		FLOATTYPE s = 123.4f;
		unsigned long j = 0;

		if(!ComputeMaxStep(i, s, j))
		{
			gLog.Print("ComputeMaxStep() failed.");

			return false;
		}

		gLog.Print("i = %ld, s = %e, j = %ld", i, s, j);

		//printf("s = %.3f, j = %ld", s, j);

		//Cycling prevention block.
		{
			static unsigned long PrevJ = -1;
			static bool WasZero = false;

			if(M_ABS(s) < M_EPSILON)//We've almost finished.
			{
				if(PrevJ == j && WasZero)
				{
					gLog.Print("##################=-PREVENTED-=##################");

					return true;
				}
			}

			PrevJ = j;
			WasZero = M_ABS(s) < M_EPSILON;
		}

		mResult += mDeltaResult * s;
		mAVector += mDeltaAVector * s;

		if(mSetC.Contains(j))
		{
			gLog.Print("Case 1");
			mSetC.Delete(j);
			mSetNC.Add(j);

			continue;
		}
		else if(mSetNC.Contains(j))
		{
			gLog.Print("Case 2");
			mSetNC.Delete(j);
			mSetC.Add(j);

			continue;
		}
		else//j must be i, implying AVector[i] = 0;
		{
			gLog.Print("Case 3");
			mSetC.Add(j);

			break;
		}
	}

	return true;
}

bool CDantzigSolver::ComputeStepDirection(unsigned long i)
{
	mDeltaResult.MakeZero();
	mDeltaResult[i] = 1.0f;

	if(mSetC.GetSize() == 0)
		return true;

	//???????????????
	CGaussSolver solver(mSetC.GetSize());

	for(unsigned long j = 0; j < mSetC.GetSize(); j++)
	{
		for(unsigned long k = 0; k < mSetC.GetSize(); k++)
			solver.mMatrix(j, k) = mAMatrix(mSetC[j], mSetC[k]);

		solver.mVector[j] = -mAMatrix(mSetC[j], i);
	}

	gLog.Print("Gauss matrix:");
	solver.mMatrix.Log();
	gLog.Print("Gauss vector:");
	solver.mVector.Log();
	gLog.Print("");

	if(!solver.Solve0())
	{
		gLog.Print("Gauss solver failed.");
		gLog.Print("Gauss matrix:");
		solver.mMatrix.Log();
		gLog.Print("");

		return false;
	}

	gLog.Print("Gauss OK:");
	solver.mResult.Log();
	gLog.Print("");

	for(unsigned long j = 0; j < mSetC.GetSize(); j++)
		mDeltaResult[mSetC[j]] = solver.mResult[j];

	return true;
}

bool CDantzigSolver::ComputeMaxStep(unsigned long i, FLOATTYPE &s, unsigned long &j)
{
	gLog.Print("----------ComputeMaxStep()----------");

	s = 1.0e36f;//A large value.
	j = -1;

	if(mDeltaAVector[i] > 0)
	{
		gLog.Print("Maxed by mDeltaAVector[i] = %e.", mDeltaAVector[i]);
		s = -mAVector[i] / mDeltaAVector[i];
		j = i;
	}

	for(unsigned long k = 0; k < mSetC.GetSize(); k++)
	{
		if(mDeltaResult[mSetC[k]] < 0)
		{
			FLOATTYPE ss = -mResult[mSetC[k]] / mDeltaResult[mSetC[k]];

			if(ss < s)
			{
				gLog.Print("Maxed by C-element: %ld = {%e\\%e}.", mSetC[k], mResult[mSetC[k]], mDeltaResult[mSetC[k]]);
				s = ss;
				j = mSetC[k];
			}
		}
	}

	for(unsigned long k = 0; k < mSetNC.GetSize(); k++)
	{
		if(mDeltaAVector[mSetNC[k]] < 0)
		{
			FLOATTYPE ss = -mAVector[mSetNC[k]] / mDeltaAVector[mSetNC[k]];

			if(ss < s)
			{
				gLog.Print("Maxed 3 by NC-element: %ld = {%e\\%e}.", mSetNC[k], mAVector[mSetNC[k]], mDeltaAVector[mSetNC[k]]);
				s = ss;
				j = mSetNC[k];
			}
		}
	}

	if(j == -1)
		return false;

	return true;
}
