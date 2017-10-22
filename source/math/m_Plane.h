#ifndef _M_PLANE_H_
#define _M_PLANE_H_

#include "all.h"

#include "m_Vector3.h"

//Класс только начат. Многие методы потребуется добавить.

/*****************************************************************************\
 * The plane is defined as ((a, b, c), d), where (a, b, c) is a unit length
 * normal vector to the plane and d is distance from the plane to the origin.
\*****************************************************************************/

class CPlane
{
public:
	union
	{
		struct
		{
			FLOATTYPE							a, b, c, d;
		};

		struct
		{
			Vec3						mNormal;		//WARNING! mNormal should be normalized for correct math!
			FLOATTYPE							mDistance;
		};
	};

	CPlane(void){}													//Do not initialize.
	CPlane(FLOATTYPE ia, FLOATTYPE ib, FLOATTYPE ic, FLOATTYPE id);
	CPlane(Vec3 normal, FLOATTYPE distance);
	CPlane(Vec3 p1, Vec3 p2, Vec3 p3);					//Make plane to go "through" 3 points.
	CPlane(const CPlane &ref);										//Copy constructor.
	//~CPlane(void){}												//Destructor. Does nothing.

	bool					operator ==(const CPlane &arg) const;	//Операция сравнения.
	bool					operator !=(const CPlane &arg) const;	//Операция сравнения.
	FLOATTYPE					GetDistanceToPoint(const Vec3 &point) const;	//Computes distance from the plane to the point.
};

#ifdef INLINE_BUILD
	#include "m_Plane.inl"
#endif

#endif


