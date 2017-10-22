#ifndef _M_INTERPOLATORS_H_
#define _M_INTERPOLATORS_H_

#include "all.h"

//Группа "сырая". адо дописать IInterpolatorable.

namespace math
{
	float LinearIP(const float &v0, const float &v1, const float &x);
	float CubicIP(const float &v0, const float &v1, const float &v2, const float &v3, const float &x);
	float BezierIP(const float &v0, const float &v1, const float &v2, const float &v3, const float &x);

	//TODO: create constructors.

	class CLinearInterpolator
	{
	public:
		float mValue[2];

		float GetPoint(const float &x) const;
	};

	class CCubicInterpolator
	{
	public:
		float mValue[4];

		float GetPoint(const float &x) const;
	};

	class CBezierInterpolator
	{
	public:
		float mValue[4];

		float GetPoint(const float &x) const;
	};
}

#endif
