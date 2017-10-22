#include <math.h>

INLINE CMatrix3::CMatrix3(void)
{
}

INLINE CMatrix3::CMatrix3(const CMatrix3 &matrix)
{
	*this = matrix;
}

INLINE CMatrix3::CMatrix3(const Vec3 &Axis, float Angle)
{
	FLOATTYPE	cf = cosf(Angle),
				sf = sinf(Angle);

	FLOATTYPE	xx = Axis.x * Axis.x,
				yy = Axis.y * Axis.y,
				zz = Axis.z * Axis.z;
	FLOATTYPE	xy = Axis.x * Axis.y,
				yz = Axis.y * Axis.z,
				zx = Axis.z * Axis.x;

	mMatrix[0][0] = xx + cf * (1.0f - xx);
	mMatrix[0][1] = xy - cf * xy - sf * Axis.z;
	mMatrix[0][2] = zx - cf * zx + sf * Axis.y;

	mMatrix[1][0] = xy - cf * xy + sf * Axis.z;
	mMatrix[1][1] = yy + cf * (1.0f - yy);
	mMatrix[1][2] = yz - cf * yz - sf * Axis.x;

	mMatrix[2][0] = zx - cf * zx - sf * Axis.y;
	mMatrix[2][1] = yz - cf * yz + sf * Axis.x;
	mMatrix[2][2] = zz + cf * (1.0f - zz);
}

/*********************************************************************************************************************\
 An zero matrix is:

  |0  0  0|
  |0  0  0|
  |0  0  0|

\*********************************************************************************************************************/
INLINE void CMatrix3::MakeZero(void)
{
	mMatrix[0][0] = 0.0;
	mMatrix[0][1] = 0.0;
	mMatrix[0][2] = 0.0;

	mMatrix[1][0] = 0.0;
	mMatrix[1][1] = 0.0;
	mMatrix[1][2] = 0.0;

	mMatrix[2][0] = 0.0;
	mMatrix[2][1] = 0.0;
	mMatrix[2][2] = 0.0;
}

/*********************************************************************************************************************\
 An identity matrix is:

  |1  0  0|
  |0  1  0|
  |0  0  1|

\*********************************************************************************************************************/
INLINE void CMatrix3::MakeIdentity(void)
{
	mMatrix[0][0] = 1.0;
	mMatrix[0][1] = 0.0;
	mMatrix[0][2] = 0.0;

	mMatrix[1][0] = 0.0;
	mMatrix[1][1] = 1.0;
	mMatrix[1][2] = 0.0;

	mMatrix[2][0] = 0.0;
	mMatrix[2][1] = 0.0;
	mMatrix[2][2] = 1.0;
}

INLINE void CMatrix3::Transpose(void)
{
	FLOATTYPE t;

	t = mMatrix[1][0];
	mMatrix[1][0] = mMatrix[0][1];
	mMatrix[0][1] = t;

	t = mMatrix[2][0];
	mMatrix[2][0] = mMatrix[0][2];
	mMatrix[0][2] = t;

	t = mMatrix[2][1];
	mMatrix[2][1] = mMatrix[1][2];
	mMatrix[1][2] = t;
}

INLINE CMatrix3 CMatrix3::GetTransposed(void) const
{
	CMatrix3 res = *this;

	res.Transpose();

	return res;
}

INLINE void CMatrix3::ToColumnMatrix4(float m[16]) const
{
	m[0]	= (float)mMatrix[0][0];
	m[1]	= (float)mMatrix[1][0];
	m[2]	= (float)mMatrix[2][0];
	m[3]	= (float)0;
	m[4]	= (float)mMatrix[0][1];
	m[5]	= (float)mMatrix[1][1];
	m[6]	= (float)mMatrix[2][1];
	m[7]	= (float)0;
	m[8]	= (float)mMatrix[0][2];
	m[9]	= (float)mMatrix[1][2];
	m[10]	= (float)mMatrix[2][2];
	m[11]	= (float)0;
	m[12]	= (float)0;
	m[13]	= (float)0;
	m[14]	= (float)0;
	m[15]	= (float)1;
}

INLINE FLOATTYPE *CMatrix3::operator [](unsigned long index)
{
	return mMatrix[index];
}

INLINE const FLOATTYPE *CMatrix3::operator[](unsigned long index) const
{
	return mMatrix[index];
}

INLINE CMatrix3 &CMatrix3::operator =(const CMatrix3 &arg)
{
	if(&arg == this)
		return *this;

	mMatrix[0][0] = arg.mMatrix[0][0];
	mMatrix[0][1] = arg.mMatrix[0][1];
	mMatrix[0][2] = arg.mMatrix[0][2];

	mMatrix[1][0] = arg.mMatrix[1][0];
	mMatrix[1][1] = arg.mMatrix[1][1];
	mMatrix[1][2] = arg.mMatrix[1][2];

	mMatrix[2][0] = arg.mMatrix[2][0];
	mMatrix[2][1] = arg.mMatrix[2][1];
	mMatrix[2][2] = arg.mMatrix[2][2];

	return *this;
}

INLINE CMatrix3 CMatrix3::operator +(const CMatrix3 &arg) const
{
	CMatrix3 result;

	for(unsigned char col = 0; col < 3; col++)
	{
		for(unsigned char row = 0; row < 3; row++)
		{
			result[row][col] = (*this)[row][col] + arg[row][col];
		}
	}

	return result;
}

INLINE CMatrix3 CMatrix3::operator -(const CMatrix3 &arg) const
{
	CMatrix3 result;

	for(unsigned char col = 0; col < 3; col++)
	{
		for(unsigned char row = 0; row < 3; row++)
		{
			result[row][col] = (*this)[row][col] - arg[row][col];
		}
	}

	return result;
}

INLINE CMatrix3 CMatrix3::operator *(FLOATTYPE arg) const
{
	CMatrix3 result;

	for(unsigned char col = 0; col < 3; col++)
	{
		for(unsigned char row = 0; row < 3; row++)
		{
			result[row][col] = (*this)[row][col] * arg;
		}
	}

	return result;
}

INLINE Vec3 CMatrix3::operator *(const Vec3 &arg) const
{
	Vec3 result;

	result.x = arg.x * mMatrix[0][0] + arg.y * mMatrix[0][1] + arg.z * mMatrix[0][2];
	result.y = arg.x * mMatrix[1][0] + arg.y * mMatrix[1][1] + arg.z * mMatrix[1][2];
	result.z = arg.x * mMatrix[2][0] + arg.y * mMatrix[2][1] + arg.z * mMatrix[2][2];

	return result;
}

INLINE CMatrix3 CMatrix3::operator *(const CMatrix3 &arg) const
{
	CMatrix3 result;

	for(unsigned char col = 0; col < 3; col++)
	{
		for(unsigned char row = 0; row < 3; row++)
		{
			result[row][col] = (*this)[row][0] * arg[0][col] +
							   (*this)[row][1] * arg[1][col] +
							   (*this)[row][2] * arg[2][col];
		}
	}

	return result;
}

INLINE void CMatrix3::operator +=(const CMatrix3 &arg)
{
	for(unsigned char col = 0; col < 3; col++)
	{
		for(unsigned char row = 0; row < 3; row++)
		{
			(*this)[row][col] += arg[row][col];
		}
	}
}

INLINE void CMatrix3::operator -=(const CMatrix3 &arg)
{
	for(unsigned char col = 0; col < 3; col++)
	{
		for(unsigned char row = 0; row < 3; row++)
		{
			(*this)[row][col] -= arg[row][col];
		}
	}
}

INLINE void CMatrix3::operator *=(FLOATTYPE arg)
{
	for(unsigned char col = 0; col < 3; col++)
	{
		for(unsigned char row = 0; row < 3; row++)
		{
			(*this)[row][col] *= arg;
		}
	}
}

INLINE void CMatrix3::operator *=(const CMatrix3 &arg)
{
	CMatrix3 result;

	for(unsigned char col = 0; col < 3; col++)
	{
		for(unsigned char row = 0; row < 3; row++)
		{
			result[row][col] = (*this)[row][0] * arg[0][col] +
							   (*this)[row][1] * arg[1][col] +
							   (*this)[row][2] * arg[2][col];
		}
	}

	*this = result;
}

INLINE CMatrix3 operator *(FLOATTYPE scalar, const CMatrix3 &matrix)
{
	return matrix * scalar;
}

INLINE Vec3 operator *(const Vec3 &vec, const CMatrix3 &mat)
{
	Vec3 result;

	result.x = vec.x * mat.mMatrix[0][0] + vec.y * mat.mMatrix[1][0] + vec.z * mat.mMatrix[2][0];
	result.y = vec.x * mat.mMatrix[0][1] + vec.y * mat.mMatrix[1][1] + vec.z * mat.mMatrix[2][1];
	result.z = vec.x * mat.mMatrix[0][2] + vec.y * mat.mMatrix[1][2] + vec.z * mat.mMatrix[2][2];

	return result;
}

INLINE FLOATTYPE CMatrix3::Det() const
{
	return mMatrix[0][0] * mMatrix[1][1] * mMatrix[2][2] +
						 mMatrix[2][0] * mMatrix[0][1] * mMatrix[1][2] +
						 mMatrix[0][2] * mMatrix[1][0] * mMatrix[2][1] -

						 mMatrix[0][2] * mMatrix[1][1] * mMatrix[2][0] -
						 mMatrix[2][2] * mMatrix[1][0] * mMatrix[0][1] -
						 mMatrix[0][0] * mMatrix[1][2] * mMatrix[2][1];
}