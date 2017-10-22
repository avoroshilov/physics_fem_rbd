#include "m_Interpolators.h"

using namespace math;

float math::LinearIP(const float &v0, const float &v1, const float &x)
{
	//return v0*(1.0f - x) + v1*x;

	return v0 + (v1 - v0) * x;
}

float math::CubicIP(const float &v0, const float &v1, const float &v2, const float &v3, const float &x)
{
	float p = (v3 - v2) - (v0 - v1);
	float q = (v0 - v1) - p;
	float r = v2 - v0;
	float s = v1;

	return p*x*x*x + q*x*x + r*x + s;
}

float math::BezierIP(const float &v0, const float &v1, const float &v2, const float &v3, const float &x)
{
	if(x <= 0)
		return v0;

	if(x >= 1)
		return v3;

	//x in power of n
	float xpow2 = x * x;
	float xpow3 = xpow2 * x;

	float ix = 1 - x;

	//1-x in power of n
	float ixpow2 = ix * ix;
	float ixpow3 = ixpow2 * ix;

	/*--------------------------------------------------------------------------------------------*\

	 r = ixpow3*v0 + 3*x*ixpow2*v1 + 3*xpow2*ix*v2 + xpow3*v3;
	     [ k0 ]	     [   k1   ]      [   k2   ]      [k3 ]

	\*--------------------------------------------------------------------------------------------*/

	/*float k0 = ixpow3;
	float k1 = 3*x*ixpow2;
	float k2 = 3*xpow2*ix;
	float k3 = xpow3;*/

	return ixpow3 * v0 + 3 * x * ixpow2 * v1 + 3 * xpow2 * ix * v2 + xpow3 * v3;
}

float math::CLinearInterpolator::GetPoint(const float &x) const
{
	return LinearIP(mValue[0], mValue[1], x);
}

float math::CCubicInterpolator::GetPoint(const float &x) const
{
	return CubicIP(mValue[0], mValue[1], mValue[2], mValue[3], x);
}

float math::CBezierInterpolator::GetPoint(const float &x) const
{
	return BezierIP(mValue[0], mValue[1], mValue[2], mValue[3], x);
}
