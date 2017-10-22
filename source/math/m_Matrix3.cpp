#include "m_Matrix3.h"

#ifndef INLINE_BUILD
	#include "m_Matrix3.inl"
#endif

/*********************************************************************************************************************\
 Algorithm uses Gram-Schmidt orthogonalization.  If 'this' matrix is M = [m0|m1|m2], then orthonormal
 output matrix is Q = [q0|q1|q2],

  q0 = m0/|m0|
  q1 = (m1-(q0*m1)q0)/|m1-(q0*m1)q0|
  q2 = (m2-(q0*m2)q0-(q1*m2)q1)/|m2-(q0*m2)q0-(q1*m2)q1|

 where |V| indicates length of vector V and A*B indicates dot product of vectors A and B.

  Думаю, что базис из столбцов матрицы получится ортонормированный. TODO: проверить.
\*********************************************************************************************************************/
void CMatrix3::Orthogonalize(void)
{
	FLOATTYPE InvLen, dot0, dot1;

	//Compute q0:
	InvLen = 1.0f / sqrt(mMatrix[0][0] * mMatrix[0][0] + mMatrix[1][0] * mMatrix[1][0] + mMatrix[2][0] * mMatrix[2][0]);

    mMatrix[0][0] *= InvLen;
    mMatrix[1][0] *= InvLen;
    mMatrix[2][0] *= InvLen;

    //Compute q1:
    dot0 = mMatrix[0][0] * mMatrix[0][1] + mMatrix[1][0] * mMatrix[1][1] + mMatrix[2][0] * mMatrix[2][1];

    mMatrix[0][1] -= dot0 * mMatrix[0][0];
    mMatrix[1][1] -= dot0 * mMatrix[1][0];
    mMatrix[2][1] -= dot0 * mMatrix[2][0];

    InvLen = 1.0f / sqrt(mMatrix[0][1] * mMatrix[0][1] + mMatrix[1][1] * mMatrix[1][1] + mMatrix[2][1] * mMatrix[2][1]);

    mMatrix[0][1] *= InvLen;
    mMatrix[1][1] *= InvLen;
    mMatrix[2][1] *= InvLen;

    //Compute q2:
    dot1 = mMatrix[0][1] * mMatrix[0][2] + mMatrix[1][1] * mMatrix[1][2] + mMatrix[2][1] * mMatrix[2][2];
	dot0 = mMatrix[0][0] * mMatrix[0][2] + mMatrix[1][0] * mMatrix[1][2] + mMatrix[2][0] * mMatrix[2][2];

    mMatrix[0][2] -= dot0 * mMatrix[0][0] + dot1 * mMatrix[0][1];
    mMatrix[1][2] -= dot0 * mMatrix[1][0] + dot1 * mMatrix[1][1];
    mMatrix[2][2] -= dot0 * mMatrix[2][0] + dot1 * mMatrix[2][1];

    InvLen = 1.0f / sqrt(mMatrix[0][2] * mMatrix[0][2] + mMatrix[1][2] * mMatrix[1][2] + mMatrix[2][2] * mMatrix[2][2]);

    mMatrix[0][2] *= InvLen;
    mMatrix[1][2] *= InvLen;
    mMatrix[2][2] *= InvLen;
}

CMatrix3 CMatrix3::GetInverse(FLOATTYPE *det)
{
	FLOATTYPE determinant = mMatrix[0][0] * mMatrix[1][1] * mMatrix[2][2] +
						 mMatrix[2][0] * mMatrix[0][1] * mMatrix[1][2] +
						 mMatrix[0][2] * mMatrix[1][0] * mMatrix[2][1] -

						 mMatrix[0][2] * mMatrix[1][1] * mMatrix[2][0] -
						 mMatrix[2][2] * mMatrix[1][0] * mMatrix[0][1] -
						 mMatrix[0][0] * mMatrix[1][2] * mMatrix[2][1];

	if(det)
		*det = determinant;

	if(M_EQUAL_AT(determinant, 0.0f, 0.00001f))
	{
		CMatrix3 res;

		res.MakeZero();

		return res;
	}

	CMatrix3 cofactors;

	cofactors[0][0] =   mMatrix[1][1] * mMatrix[2][2] - mMatrix[2][1] * mMatrix[1][2];
	cofactors[0][1] = -(mMatrix[1][0] * mMatrix[2][2] - mMatrix[2][0] * mMatrix[1][2]);
	cofactors[0][2] =   mMatrix[1][0] * mMatrix[2][1] - mMatrix[2][0] * mMatrix[1][1];

	cofactors[1][0] = -(mMatrix[0][1] * mMatrix[2][2] - mMatrix[2][1] * mMatrix[0][2]);
	cofactors[1][1] =   mMatrix[0][0] * mMatrix[2][2] - mMatrix[2][0] * mMatrix[0][2];
	cofactors[1][2] = -(mMatrix[0][0] * mMatrix[2][1] - mMatrix[2][0] * mMatrix[0][1]);

	cofactors[2][0] =   mMatrix[0][1] * mMatrix[1][2] - mMatrix[1][1] * mMatrix[0][2];
	cofactors[2][1] = -(mMatrix[0][0] * mMatrix[1][2] - mMatrix[1][0] * mMatrix[0][2]);
	cofactors[2][2] =   mMatrix[0][0] * mMatrix[1][1] - mMatrix[1][0] * mMatrix[0][1];

	cofactors.Transpose();

	return (1.0f / determinant) * cofactors;
}

CMatrix3 CMatrix3::Sqrt(unsigned int niter) const
{
	CMatrix3 root = *this, invroot, temp;
	invroot.MakeIdentity();

	for(unsigned int i = 0; i < niter; ++i)
	{
		temp = root;
		root = 0.5f  * (root + invroot.GetInverse());
		invroot = 0.5f * (invroot + temp.GetInverse());

		// Newton method (unstable)
//		root = 0.5f  * (root + root.GetInverse() * *this);
	}

	return root;
}

void CMatrix3::PolarDecompose(CMatrix3 &p, CMatrix3 &u, unsigned int niter) const
{
	p = (GetTransposed() * *this).Sqrt(niter);
	u = *this * p.GetInverse();
}
