#include <math.h>

INLINE CMatrix4::CMatrix4(void)
{
	this->MakeIdentity();

	return;
}

INLINE CMatrix4::CMatrix4(const CMatrix4 &matrix)
{
	*this = matrix;

	return;
}

INLINE CMatrix4::CMatrix4(const CMatrix3 &rotate, const Vec3 &translate)
{
	mMatrix[0][0] = rotate.mMatrix[0][0];
	mMatrix[0][1] = rotate.mMatrix[0][1];
	mMatrix[0][2] = rotate.mMatrix[0][2];

	mMatrix[1][0] = rotate.mMatrix[1][0];
	mMatrix[1][1] = rotate.mMatrix[1][1];
	mMatrix[1][2] = rotate.mMatrix[1][2];

	mMatrix[2][0] = rotate.mMatrix[2][0];
	mMatrix[2][1] = rotate.mMatrix[2][1];
	mMatrix[2][2] = rotate.mMatrix[2][2];

	mMatrix[0][3] = translate[0];
	mMatrix[1][3] = translate[1];
	mMatrix[2][3] = translate[2];

	mMatrix[3][0] = 0.0;
	mMatrix[3][1] = 0.0;
	mMatrix[3][2] = 0.0;
	mMatrix[3][3] = 1.0;

	return;
}

//=====================================================================================================================
// A zero matrix is:
//
//  |0  0  0  0|
//  |0  0  0  0|
//  |0  0  0  0|
//  |0  0  0  0|
//
//=====================================================================================================================
INLINE void CMatrix4::MakeZero(void)
{
	for(unsigned long row = 0; row < 4; row++)
	{
		for(unsigned long col = 0; col < 4; col++)
		{
			mMatrix[row][col] = 0.0;
		}
	}

	return;
}

//=====================================================================================================================
// An identity matrix is:
//
//  |1  0  0  0|
//  |0  1  0  0|
//  |0  0  1  0|
//  |0  0  0  1|
//
//=====================================================================================================================
INLINE void CMatrix4::MakeIdentity(void)
{
	for(unsigned long row = 0; row < 4; row++)
	{
		for(unsigned long col = 0; col < 4; col++)
		{
			mMatrix[row][col] = row == col ? 1.0f : 0.0f;
		}
	}

	return;
}

//=====================================================================================================================

INLINE void CMatrix4::MakeZooming(const Vec3 &origin, const FLOATTYPE &zoom)
{
	FLOATTYPE z = 1.0f - zoom;

	mMatrix[0][0] = zoom;
	mMatrix[0][1] = 0;
	mMatrix[0][2] = 0;

	mMatrix[1][0] = 0;
	mMatrix[1][1] = zoom;
	mMatrix[1][2] = 0;

	mMatrix[2][0] = 0;
	mMatrix[2][1] = 0;
	mMatrix[2][2] = zoom;

	mMatrix[3][0] = 0;
	mMatrix[3][1] = 0;
	mMatrix[3][2] = 0;

	mMatrix[0][3] = z * origin.x;
	mMatrix[1][3] = z * origin.y;
	mMatrix[2][3] = z * origin.z;
	mMatrix[3][3] = 1;

	return;
}

//=====================================================================================================================

INLINE FLOATTYPE *CMatrix4::operator [](unsigned long index)
{
	return mMatrix[index];
}

//=====================================================================================================================

INLINE const FLOATTYPE *CMatrix4::operator[](unsigned long index) const
{
	return mMatrix[index];
}

//=====================================================================================================================

INLINE CMatrix4 &CMatrix4::operator =(const CMatrix4 &arg)
{
	if(&arg == this)
		return *this;

	for(unsigned long row = 0; row < 4; row++)
	{
		for(unsigned long col = 0; col < 4; col++)
		{
			mMatrix[row][col] = arg.mMatrix[row][col];
		}
	}

	return *this;
}

//=====================================================================================================================

INLINE CMatrix4 CMatrix4::operator *(FLOATTYPE arg) const
{
	CMatrix4 result;

	for(unsigned long row = 0; row < 4; row++)
	{
		for(unsigned long col = 0; col < 4; col++)
		{
			result.mMatrix[row][col] = mMatrix[row][col] * arg;
		}
	}

	return result;
}

//=====================================================================================================================

INLINE Vec3 CMatrix4::operator *(const Vec3 &arg) const
{
	Vec3 result;

	result.x = arg.x * mMatrix[0][0] + arg.y * mMatrix[0][1] + arg.z * mMatrix[0][2] + mMatrix[0][3];
	result.y = arg.x * mMatrix[1][0] + arg.y * mMatrix[1][1] + arg.z * mMatrix[1][2] + mMatrix[1][3];
	result.z = arg.x * mMatrix[2][0] + arg.y * mMatrix[2][1] + arg.z * mMatrix[2][2] + mMatrix[2][3];

	return result;
}

//=====================================================================================================================

INLINE CMatrix4 CMatrix4::operator *(const CMatrix4 &arg) const
{
	//64 muls + 48 adds

	CMatrix4 result;

	for(unsigned char row = 0; row < 4; row++)
	{
		for(unsigned char col = 0; col < 4; col++)
		{
			result.mMatrix[row][col] =
				mMatrix[row][0] * arg.mMatrix[0][col] +
				mMatrix[row][1] * arg.mMatrix[1][col] +
				mMatrix[row][2] * arg.mMatrix[2][col] +
				mMatrix[row][3] * arg.mMatrix[3][col];
		}
	}

	return result;
}

//=====================================================================================================================

INLINE CMatrix4 operator *(FLOATTYPE scalar, const CMatrix4 &matrix)
{
	return matrix * scalar;
}

//=====================================================================================================================

INLINE Vec3 CMatrix4::RotateOnly(const Vec3 &arg) const
{
	Vec3 result;

	result.x = arg.x * mMatrix[0][0] + arg.y * mMatrix[0][1] + arg.z * mMatrix[0][2];
	result.y = arg.x * mMatrix[1][0] + arg.y * mMatrix[1][1] + arg.z * mMatrix[1][2];
	result.z = arg.x * mMatrix[2][0] + arg.y * mMatrix[2][1] + arg.z * mMatrix[2][2];

	return result;
}

//=====================================================================================================================

INLINE CMatrix4 &CMatrix4::Add(const CMatrix4 &arg)
{
	for(unsigned char row = 0; row < 4; row++)
	{
		for(unsigned char col = 0; col < 4; col++)
		{
			mMatrix[row][col] += arg.mMatrix[row][col];
		}
	}

	return *this;
}

//=====================================================================================================================

INLINE void CMatrix4::MulTranslationByIt(const Vec3 &t)
{
	/*

		|1 0 0 x|   |a b c d|   |a+nx b+px c+qx d+sx|
		|0 1 0 y| * |e f g h| = |e+ny f+py g+qy h+sy|
		|0 0 1 z|   |i j k m|   |i+nz j+pz k+qz m+sz|
		|0 0 0 1|   |n p q s|   |n    p    q    s   |

		12 muls
		12 adds (+4 adds)

	*/

	CMatrix4 m;

	m.mMatrix[0][0] = mMatrix[3][0] * t.x;
	m.mMatrix[0][1] = mMatrix[3][1] * t.x;
	m.mMatrix[0][2] = mMatrix[3][2] * t.x;
	m.mMatrix[0][3] = mMatrix[3][3] * t.x;

	m.mMatrix[1][0] = mMatrix[3][0] * t.y;
	m.mMatrix[1][1] = mMatrix[3][1] * t.y;
	m.mMatrix[1][2] = mMatrix[3][2] * t.y;
	m.mMatrix[1][3] = mMatrix[3][3] * t.y;

	m.mMatrix[2][0] = mMatrix[3][0] * t.z;
	m.mMatrix[2][1] = mMatrix[3][1] * t.z;
	m.mMatrix[2][2] = mMatrix[3][2] * t.z;
	m.mMatrix[2][3] = mMatrix[3][3] * t.z;

	m.mMatrix[3][0] = 0;
	m.mMatrix[3][1] = 0;
	m.mMatrix[3][2] = 0;
	m.mMatrix[3][3] = 0;

	this->Add(m);
	
	return;
}

//=====================================================================================================================

INLINE void CMatrix4::MulItByTranslation(const Vec3 &t)
{
	/*

		|a b c d|   |1 0 0 x|   |a b c ax+by+cz+d|
		|e f g h| * |0 1 0 y| = |e f g ex+fy+gz+h|
		|i j k m|   |0 0 1 z|   |i j k ix+jy+kz+m|
		|n p q s|   |0 0 0 1|   |n p q nx+py+qz+s|

		12 muls
		16 adds (+12 adds)

	*/

	FLOATTYPE v[4];

	v[0] = mMatrix[0][0] + t.x + mMatrix[0][1] * t.y + mMatrix[0][2] * t.z;
	v[1] = mMatrix[1][0] + t.x + mMatrix[1][1] * t.y + mMatrix[1][2] * t.z;
	v[2] = mMatrix[2][0] + t.x + mMatrix[2][1] * t.y + mMatrix[2][2] * t.z;
	v[3] = mMatrix[3][0] + t.x + mMatrix[3][1] * t.y + mMatrix[3][2] * t.z;

	mMatrix[0][3] += v[0];
	mMatrix[1][3] += v[1];
	mMatrix[2][3] += v[2];
	mMatrix[3][3] += v[3];

	return;
}
