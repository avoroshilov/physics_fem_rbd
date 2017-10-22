#ifndef _M_MATH_FUNCTIONS_H_
#define _M_MATH_FUNCTIONS_H_

#include "all.h"
#include <vector>

#include "m_Vector3.h"

namespace math
{
	// Computes a normal to arbitrary 3D polygon
	Vec3									NewellNormalize(const std::vector<Vec3> &perimeter);

	// Does the face look at the point?
	bool									LooksAt(Vec3 LookPoint, Vec3 p1, Vec3 p2, Vec3 p3);

	// Does the face look at the point? A variant for optimizing.
	bool									LooksAtB(Vec3 LookPoint, Vec3 p1, Vec3 p2, Vec3 p3);

	// Distance to triangle's surface.
	FLOATTYPE								Distance(Vec3 P, Vec3 p1, Vec3 p2, Vec3 p3);

	FLOATTYPE								Factorial(unsigned long n);

	FLOATTYPE								CMN(unsigned long m, unsigned long n);

	// Computes derivative over time of distance between two lines.
	// Lines are defined by two segments.
	// Considering the segments as edges of two moving objects.
	//
	// r1, r2, r3, r4 - four points defining two segments (two tail-head pairs)
	// opos1, opos2   - objects' positions
	// vel1, vel2     - objects' linear velocities
	// omega1, omega2 - objects' angular velocities
	//
	// pDistSqrDOT and pDist are optional.
	FLOATTYPE								GetLinesDistanceDOT(const Vec3 &r1, const Vec3 &r2,
											const Vec3 &r3, const Vec3 &r4,
											const Vec3 &opos1, const Vec3 &opos2,
											const Vec3 &vel1, const Vec3 &vel2,
											const Vec3 &omega1, const Vec3 &omega2,
											FLOATTYPE *pDistSqrDOT = NULL, FLOATTYPE *pDist = NULL);

	// Computes derivative over time of distance between a vertex and a face.
	// Considering the vertex and the face as features on two moving objects.
	FLOATTYPE								GetVFDistanceDOT(const Vec3 &r,
											const Vec3 &v,//a vertex from the face
											const Vec3 &vv,//transformed (v + normal)
											const Vec3 &opos1, const Vec3 &opos2,
											const Vec3 &vel1, const Vec3 &vel2,
											const Vec3 &omega1, const Vec3 &omega2);
	
	//returns a vector, orthogonal to the given (one of many). If the given vector is
	//degenerate, degenerate vector is returned.
	Vec3									ChooseRandomOrthogonal(const Vec3 &arg);
	

	// Fast absolute value (sign bit zeroed)
	__forceinline float						fastabs(float f)
	{
		int i = ((*(int *)&f)&0x7fffffff);
		return (*(float *)&i);
	}
}



#endif
