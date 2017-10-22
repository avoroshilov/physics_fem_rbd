#include <math.h>

#include "m_Const.h"

#include "helpers/debug.h"

INLINE Vec3::Vec3(FLOATTYPE ix, FLOATTYPE iy, FLOATTYPE iz)
{
	x = ix; y = iy; z = iz;
}

INLINE void Vec3::Assign(FLOATTYPE ix, FLOATTYPE iy, FLOATTYPE iz)
{
	x = ix;
	y = iy;
	z = iz;
}

INLINE void Vec3::Set(FLOATTYPE ix, FLOATTYPE iy, FLOATTYPE iz)
{
	x = ix;
	y = iy;
	z = iz;
}

INLINE FLOATTYPE Vec3::GetLength(void) const
{
	return sqrt(x*x + y*y + z*z);
}

INLINE FLOATTYPE Vec3::GetSquareLength(void) const
{
	return x*x + y*y + z*z;
}

INLINE void Vec3::MakeInverse(void)
{
	x = -x;
	y = -y;
	z = -z;
}

INLINE FLOATTYPE Vec3::Dot(const Vec3 &arg) const
{
	return x*arg.x + y*arg.y + z*arg.z;
}

INLINE Vec3 Vec3::Cross(const Vec3 &arg) const
{
	return Vec3(y*arg.z - z*arg.y, z*arg.x - x*arg.z, x*arg.y - y*arg.x);
}

INLINE FLOATTYPE Vec3::Triple(const Vec3 &arg1, const Vec3 &arg2) const
{
	// ax*(by*cz - cy*bz) + ay*(cx*bz - bx*cz) + az*(bx*cy - cx*by)
	return x * (arg1.y * arg2.z - arg2.y * arg1.z) + y * (arg2.x * arg1.z - arg1.x * arg2.z) + z * (arg1.x * arg2.y - arg2.x * arg1.y);
}

INLINE Vec3 Vec3::ProjectOnto(const Vec3 &arg) const
{
	return arg * this->Dot(arg) / arg.GetSquareLength();
}

INLINE FLOATTYPE &Vec3::operator [](unsigned long index)
{
	return v[index];
}

INLINE const FLOATTYPE &Vec3::operator[](unsigned long index) const
{
	return v[index];
}

INLINE const Vec3 &Vec3::operator =(const Vec3 &arg)
{
	if(this == &arg)
	{
		return *this;
	}

	x = arg.x;
	y = arg.y;
	z = arg.z;

	return *this;
}

INLINE bool Vec3::operator ==(const Vec3 &arg) const
{
	return (M_EQUAL(x, arg.x) && M_EQUAL(y, arg.y) && M_EQUAL(z, arg.z));
}

INLINE bool Vec3::operator !=(const Vec3 &arg) const
{
	return !(*this == arg);
}

INLINE bool Vec3::IsEqualAt(const Vec3 &arg, const FLOATTYPE &precision) const
{
	return (M_EQUAL_AT(x, arg.x, precision) && M_EQUAL_AT(y, arg.y, precision) && M_EQUAL_AT(z, arg.z, precision));
}

INLINE bool Vec3::operator <(const Vec3 &arg) const
{
	if(x == arg.x && y == arg.y)
	{
		return z < arg.z;
	}

	if(x == arg.x)
	{
		return y < arg.y;
	}

	return x < arg.x;
}

INLINE Vec3 Vec3::operator -(void) const//monadic -
{
	return Vec3(-x, -y, -z);
}

INLINE Vec3 Vec3::operator +(const Vec3 &arg) const
{
	return Vec3(x + arg.x, y + arg.y, z + arg.z);
}

INLINE Vec3 Vec3::operator -(const Vec3 &arg) const//binary -
{
	return Vec3(x - arg.x, y - arg.y, z - arg.z);
}

INLINE Vec3 Vec3::operator *(const FLOATTYPE &arg) const
{
	return Vec3(x * arg, y * arg, z * arg);
}

INLINE Vec3 Vec3::operator /(const FLOATTYPE &arg) const
{
	return Vec3(x / arg, y / arg, z / arg);
}

INLINE Vec3 &Vec3::operator +=(const Vec3 &arg)
{
	x += arg.x;
	y += arg.y;
	z += arg.z;

	return *this;
}

INLINE Vec3 &Vec3::operator -=(const Vec3 &arg)
{
	x -= arg.x;
	y -= arg.y;
	z -= arg.z;

	return *this;
}

INLINE Vec3 &Vec3::operator *=(const FLOATTYPE &arg)
{
	x *= arg;
	y *= arg;
	z *= arg;

	return *this;
}

INLINE Vec3 &Vec3::operator /=(const FLOATTYPE &arg)
{
	x /= arg;
	y /= arg;
	z /= arg;

	return *this;
}

INLINE Vec3 operator *(const FLOATTYPE &scalar, const Vec3 &vector)
{
	return vector * scalar;
}
