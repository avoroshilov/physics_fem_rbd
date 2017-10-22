#include "helpers/debug.h"
#include "m_Const.h"
#include <math.h>

INLINE CPlane::CPlane(FLOATTYPE ia, FLOATTYPE ib, FLOATTYPE ic, FLOATTYPE id)
{
	assert(M_EQUAL(Vec3(ia, ib, ic).GetSquareLength(), 1.0));//a bug fixed

	a = ia;
	b = ib;
	c = ic;
	d = id;
}

INLINE CPlane::CPlane(Vec3 normal, FLOATTYPE distance)
{
	assert(M_ABS(normal.GetSquareLength() - 1.0) < M_EPSILON);

	mNormal = normal;
	mDistance = distance;
}

INLINE CPlane::CPlane(Vec3 p1, Vec3 p2, Vec3 p3)
{
	Vec3 v1, v2;

	v1 = p2 - p1;
	v2 = p3 - p1;

	a = v1.y * v2.z - v2.y * v1.z;
	b = v1.z * v2.x - v2.z * v1.x;
	c = v1.x * v2.y - v2.x * v1.y;

	mNormal.Normalize();

	d = -mNormal.Dot(p1);
}

INLINE CPlane::CPlane(const CPlane &ref)
{
	mNormal = ref.mNormal;
	mDistance = ref.mDistance;
}

INLINE bool CPlane::operator ==(const CPlane &arg) const
{
	return (M_EQUAL_AT(a, arg.a, 0.01 * M_EPSILON)) && (M_EQUAL_AT(b, arg.b, 0.01 * M_EPSILON)) && (M_EQUAL_AT(c, arg.c, 0.01 * M_EPSILON)) && (M_EQUAL_AT(d, arg.d, 0.01 * M_EPSILON) );
}

INLINE bool	CPlane::operator !=(const CPlane &arg) const
{
	return !(*this == arg);
}

INLINE FLOATTYPE CPlane::GetDistanceToPoint(const Vec3 &point) const
{
	return mNormal.Dot(point) + mDistance;
}
