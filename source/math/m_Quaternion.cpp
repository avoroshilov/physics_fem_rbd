#include "m_Quaternion.h"

const unsigned char CQuaternion::mNext[3] = {1, 2, 0};

#ifndef INLINE_BUILD
	#include "m_Quaternion.inl"
#endif

void CQuaternion::MakeInverse(void)
{
	CQuaternion inverse;

	FLOATTYPE norm = 0.0;

	for(unsigned char i = 0; i < 4; i++)
		norm += mCoords[i] * mCoords[i];

	if(norm > M_EPSILON)
	{
		FLOATTYPE InvNorm = 1.0f / norm;
		inverse.w =  w * InvNorm;
		inverse.x = -x * InvNorm;
		inverse.y = -y * InvNorm;
		inverse.z = -z * InvNorm;
	}
	else
	{
		message("Quaternion dead.");
		//Return an invalid result to flag the error
		for(unsigned char i = 0; i < 4; i++)
			inverse.mCoords[i] = 0.0;
	}

	*this = inverse;
}

void CQuaternion::FromMatrix3(const CMatrix3 &matrix)
{
/*********************************************************************************************************************\
  Algorithm from Ken Shoemake's article in 1987 SIGGRAPH course notesarticle "Quaternion Calculus and Fast Animation".
\*********************************************************************************************************************/
	FLOATTYPE trace = matrix[0][0] + matrix[1][1] + matrix[2][2];
	FLOATTYPE root;

	if(trace > 0.0)
	{
		//|w| > 1/2, may as well choose w > 1/2
		root = sqrt(trace + 1.0f);//2w
		w = 0.5f * root;
		root = 0.5f / root;//1/(4w)

		x = (matrix[2][1] - matrix[1][2]) * root;
		y = (matrix[0][2] - matrix[2][0]) * root;
		z = (matrix[1][0] - matrix[0][1]) * root;
	}
	else
	{
		//|w| <= 1/2
		int i = 0;

		if(matrix[1][1] > matrix[0][0])
			i = 1;

		if(matrix[2][2] > matrix[i][i])
			i = 2;

		int j = mNext[i];
		int k = mNext[j];

		root = sqrt(matrix[i][i] - matrix[j][j] - matrix[k][k] + 1.0f);

		//FLOATTYPE *q[3] = { &x, &y, &z };
		//*q[i] = 0.5 * root;
		mCoords[i + 1] = 0.5f * root;//EDITED 2006.11.09, 15.41, "+ 1' added.

		root = 0.5f / root;
		w = (matrix[k][j] - matrix[j][k]) * root;
		/* *q[j] */ mCoords[j] = (matrix[j][i] + matrix[i][j]) * root;
		/* *q[k] */ mCoords[k] = (matrix[k][i] + matrix[i][k]) * root;
	}
}

void CQuaternion::ToMatrix3(CMatrix3 &matrix) const
{
/*********************************************************************************************************************\
   The matrix corresponding to the quaternion [w, x, y, z] will be:
	 |1     - 2*y*y - 2*z*z      2*x*y - 2*w*z              2*x*z + 2*w*y        |
     |2*x*y + 2*w*z              1     - 2*x*x - 2*z*z      2*y*z - 2*w*x        |
     |2*x*z - 2*w*y              2*y*z + 2*w*x              1     - 2*x*x - 2*y*y|
\*********************************************************************************************************************/


	//TODO: REWRITE!!!!!!!!!!!!!!


	//[dx] = [FLOATTYPE x]
	FLOATTYPE dx  = 2.0f * x;
	FLOATTYPE dy  = 2.0f * y;
	FLOATTYPE dz  = 2.0f * z;

	//[dxx] = [FLOATTYPE square x] = 2 * x * x
	//The same about others...
	FLOATTYPE dwx = dx * w;
	FLOATTYPE dwy = dy * w;
	FLOATTYPE dwz = dz * w;

	FLOATTYPE dxy = dx * y;
	FLOATTYPE dxz = dx * z;
	FLOATTYPE dyz = dy * z;

	FLOATTYPE dxx = dx * x;
	FLOATTYPE dyy = dy * y;
	FLOATTYPE dzz = dz * z;

	matrix[0][0] = 1.0f - dyy - dzz;
	matrix[0][1] = dxy  - dwz;
	matrix[0][2] = dxz  + dwy;

	matrix[1][0] = dxy  + dwz;
	matrix[1][1] = 1.0f - dxx - dzz;
	matrix[1][2] = dyz  - dwx;

	matrix[2][0] = dxz  - dwy;
	matrix[2][1] = dyz  + dwx;
	matrix[2][2] = 1.0f - dxx - dyy;
}
