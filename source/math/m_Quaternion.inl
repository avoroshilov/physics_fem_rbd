#include <math.h>
#include <memory.h>

#include "helpers/debug.h"
#include "helpers/log.h"

#include "m_Const.h"

INLINE CQuaternion::CQuaternion(FLOATTYPE iw, FLOATTYPE ix, FLOATTYPE iy, FLOATTYPE iz)
{
	w = iw;
	x = ix;
	y = iy;
	z = iz;
}

INLINE CQuaternion::CQuaternion(const CQuaternion &quat)
{
	w = quat.w;
	x = quat.x;
	y = quat.y;
	z = quat.z;
}

INLINE CQuaternion::CQuaternion(const CMatrix3 &matrix)
{
	this->FromMatrix3(matrix);
}

INLINE CQuaternion::CQuaternion(const Vec3 &axis, FLOATTYPE angle)
{
	this->FromAxisAngle(axis, angle);
}

INLINE void CQuaternion::Assign(FLOATTYPE iw, FLOATTYPE ix, FLOATTYPE iy, FLOATTYPE iz)
{
	w = iw;
	x = ix;
	y = iy;
	z = iz;
}

INLINE void CQuaternion::MakeZero(void)
{
	w = 0.0;
	x = 0.0;
	y = 0.0;
	z = 0.0;
}

INLINE void CQuaternion::MakeIdentity(void)
{
	w = 1.0;
	x = 0.0;
	y = 0.0;
	z = 0.0;
}

INLINE void CQuaternion::Normalize(void)
{
	FLOATTYPE magnitude = sqrt(w*w + x*x + y*y + z*z);

	assert(magnitude > M_EPSILON);

	w /= magnitude;
	x /= magnitude;
	y /= magnitude;
	z /= magnitude;
}

INLINE void CQuaternion::Conjugate(void)
{
	//Assert:  'this' is unit length

	//x = -x;
	//y = -y;
	//z = -z;

	//or

	w = -w;
}

INLINE FLOATTYPE CQuaternion::Dot(const CQuaternion &arg) const
{
	return w * arg.w + x * arg.x + y * arg.y + z * arg.z;
}

INLINE FLOATTYPE &CQuaternion::operator [](int index)
{
	assert(index >= 0 && index <= 3);

	return mCoords[index];
}

INLINE const FLOATTYPE &CQuaternion::operator [](int index) const
{
	assert(index >= 0 && index <= 3);

	return mCoords[index];
}

INLINE CQuaternion &CQuaternion::operator =(const CQuaternion &arg)
{
	if(&arg == this)
		return *this;

	memcpy(this->mCoords, arg.mCoords, 4 * sizeof(mCoords[0]));

	return *this;
}

INLINE bool CQuaternion::operator ==(const CQuaternion &arg) const
{
	for(unsigned char i = 0; i < 4; i++)
	{
		if(mCoords[i] != arg.mCoords[i])
			return false;
	}

	return true;
}

INLINE bool CQuaternion::operator !=(const CQuaternion &arg) const
{
	return !operator ==(arg);
}

INLINE CQuaternion CQuaternion::operator -(void) const
{
	CQuaternion res;

	for(unsigned char i = 0; i < 4; i++)
		res.mCoords[i] = -mCoords[i];

	return res;
}

INLINE CQuaternion CQuaternion::operator +(const CQuaternion &arg) const
{
	CQuaternion sum;

	for(unsigned char i = 0; i < 4; i++)
		sum.mCoords[i] = mCoords[i] + arg.mCoords[i];

	return sum;
}

INLINE CQuaternion CQuaternion::operator -(const CQuaternion &arg) const
{
	CQuaternion diff;

	for(unsigned char i = 0; i < 4; i++)
		diff.mCoords[i] = mCoords[i] - arg.mCoords[i];

	return diff;
}

INLINE CQuaternion CQuaternion::operator *(const CQuaternion &arg) const
{
	//NOTE:  Multiplication is not generally commutative, so in most
	//cases p*q != q*p.

	CQuaternion kProd;

	kProd.mCoords[0] =
		mCoords[0] * arg.mCoords[0] -
		mCoords[1] * arg.mCoords[1] -
		mCoords[2] * arg.mCoords[2] -
		mCoords[3] * arg.mCoords[3];

	kProd.mCoords[1] =
		mCoords[0] * arg.mCoords[1] +
		mCoords[1] * arg.mCoords[0] +
		mCoords[2] * arg.mCoords[3] -
		mCoords[3] * arg.mCoords[2];

	kProd.mCoords[2] =
		mCoords[0] * arg.mCoords[2] +
		mCoords[2] * arg.mCoords[0] +
		mCoords[3] * arg.mCoords[1] -
		mCoords[1] * arg.mCoords[3];

	kProd.mCoords[3] =
		mCoords[0] * arg.mCoords[3] +
		mCoords[3] * arg.mCoords[0] +
		mCoords[1] * arg.mCoords[2] -
		mCoords[2] * arg.mCoords[1];

	return kProd;
}

INLINE CQuaternion CQuaternion::operator *(FLOATTYPE scalar) const
{
	CQuaternion kProd;

	for(unsigned char i = 0; i < 4; i++)
		kProd.mCoords[i] = scalar * mCoords[i];

	return kProd;
}

INLINE CQuaternion CQuaternion::operator /(FLOATTYPE scalar) const
{
	CQuaternion res;

	assert(M_ABS(scalar) > M_EPSILON && "Attempt to divide a quaternion by a small value.");

	if(scalar != 0.0)
	{
		FLOATTYPE InvScalar = 1.0f / scalar;
		for(unsigned char i = 0; i < 4; i++)
			res.mCoords[i] *= InvScalar;
	}

	return res;
}

INLINE Vec3 CQuaternion::operator *(const Vec3 &vector) const
{
	//Given a vector u =(x0, y0, z0) and a unit length quaternion
	//q = <w, x, y, z>, the vector v =(x1, y1, z1) which represents the
	//rotation of u by q is v = q*u*q^{-1} where * indicates quaternion
	//multiplication and where u is treated as the quaternion <0, x0, y0, z0>.
	//Note that q^{-1} = <w, -x, -y, -z>, so no real work is required to
	//invert q.  Now
	//
	//  q*u*q^{-1} = q*<0, x0, y0, z0>*q^{-1}
	//	= q*(x0*i+y0*j+z0*k)*q^{-1}
	//	= x0*(q*i*q^{-1})+y0*(q*j*q^{-1})+z0*(q*k*q^{-1})
	//
	//As 3-vectors, q*i*q^{-1}, q*j*q^{-1}, and q*k*q^{-1} are the columns
	//of the rotation matrix computed in CQuaternion::ToRotationMatrix.
	//The vector v is obtained as the product of that rotation matrix with
	//vector u.  As such, the quaternion representation of a rotation
	//matrix requires less space than the matrix and more time to compute
	//the rotated vector. Typical space-time tradeoff...

	CMatrix3 matrix;

	this->ToMatrix3(matrix);

	return matrix * vector;
}

INLINE CQuaternion &CQuaternion::operator +=(const CQuaternion &arg)
{
	for(unsigned char i = 0; i < 4; i++)
		mCoords[i] += arg.mCoords[i];

	return *this;
}

INLINE CQuaternion &CQuaternion::operator -=(const CQuaternion &arg)
{
	for(unsigned char i = 0; i < 4; i++)
		mCoords[i] -= arg.mCoords[i];

	return *this;
}

INLINE CQuaternion &CQuaternion::operator *=(FLOATTYPE scalar)
{
	for(unsigned char i = 0; i < 4; i++)
		mCoords[i] *= scalar;

	return *this;
}

INLINE CQuaternion &CQuaternion::operator /=(FLOATTYPE scalar)
{
	if(scalar != 0.0)
	{
		FLOATTYPE InvScalar = 1.0f / scalar;
		for(unsigned char i = 0; i < 4; i++)
			mCoords[i] *= InvScalar;
	}
	else
	{
		FLOATTYPE InvScalar = 1.0f / scalar;
		for(unsigned char i = 0; i < 4; i++)
			mCoords[i] *= InvScalar;
	}

	return *this;
}

INLINE void CQuaternion::FromAxisAngle(const Vec3 &axis, FLOATTYPE angle)
{
	//assert:  axis[] is unit length
	//
	//The quaternion representing the rotation is
	//  q = cos(A/2) + sin(A/2) * (x*i + y*j + z*k)

	FLOATTYPE halfangle = 0.5f * angle;
	FLOATTYPE ssin = sin(halfangle);
	w = cos(halfangle);
	x = ssin * axis[0];
	y = ssin * axis[1];
	z = ssin * axis[2];
}

INLINE void CQuaternion::ToAxisAngle(Vec3 &axis, FLOATTYPE &angle) const
{
	//The quaternion representing the rotation is
	//  q = cos(A/2)+sin(A/2)*(x*i+y*j+z*k)

	FLOATTYPE SqrLength = mCoords[1] * mCoords[1] + mCoords[2] * mCoords[2] + mCoords[3] * mCoords[3];

	if(SqrLength > 0.0)
	{
		angle = 2.0f * acos(mCoords[0]);
		FLOATTYPE InvLength = 1 / sqrt(SqrLength);
		axis[0] = mCoords[1] * InvLength;
		axis[1] = mCoords[2] * InvLength;
		axis[2] = mCoords[3] * InvLength;
	}
	else
	{
		//angle is 0(mod 2*pi), so any axis will do
		angle   = 0.0;
		axis[0] = 1.0;
		axis[1] = 0.0;
		axis[2] = 0.0;
	}
}

/* Dead ZET's version

INLINE CQuaternion CQuaternion::Slerp(FLOATTYPE t, const CQuaternion &p, const CQuaternion &q)
{
	FLOATTYPE ccos = p.Dot(q);
	FLOATTYPE angle = acos(ccos);

	if(M_ABS(angle) < M_EPSILON)
		return p;

	FLOATTYPE ssin = sin(angle);
	FLOATTYPE InvSin = 1.0 / ssin;
	FLOATTYPE Coeff0 = sin((1.0 - t) * angle) * InvSin;
	FLOATTYPE Coeff1 = sin(t * angle) * InvSin;

	return p * Coeff0 + q * Coeff1;
}
*/

INLINE CQuaternion CQuaternion::Slerp(FLOATTYPE t, const CQuaternion &p, const CQuaternion &q)
{
	FLOATTYPE ccos = p.Dot(q);
	FLOATTYPE Coeff0 = 1.0, Coeff1 = 1.0;
	
	if(ccos < 0)
	{
		ccos = -ccos;
		Coeff1 = -1.0;
	}

	FLOATTYPE angle = acos(ccos);

	if(M_ABS(angle) < M_EPSILON)
		return p;

	FLOATTYPE ssin = sin(angle);
	FLOATTYPE InvSin = 1.0f / ssin;
	
	Coeff0 *= sin((1.0f - t) * angle) * InvSin;
	Coeff1 *= sin(t * angle) * InvSin;

	return p * Coeff0 + q * Coeff1;
}

INLINE CQuaternion CQuaternion::Squad(FLOATTYPE t, const CQuaternion &q0, const CQuaternion &a0, const CQuaternion &a1, const CQuaternion &q1)
{
	FLOATTYPE SlerpT = 2.0f * t * (1.0f - t);
	CQuaternion kSlerpP = Slerp(t, q0, q1);
	CQuaternion kSlerpQ = Slerp(t, a0, a1);

	return Slerp(SlerpT, kSlerpP, kSlerpQ);
}

INLINE CQuaternion operator *(FLOATTYPE scalar, const CQuaternion &quat)
{
	return quat * scalar;
}
