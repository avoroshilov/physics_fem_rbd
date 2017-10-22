#ifndef _GJK_H_
#define _GJK_H_

#include "math/m_Vector3.h"
#include "math/m_Matrix3.h"

class GJKObj
{
public:
	virtual ~GJKObj() 
	{
	};
	
	//direction must have a unit length
	virtual Vec3 GetSupportPoint(const Vec3 &dir) = 0;
	virtual Vec3 GetLastPoint() = 0;

};

float GJK(GJKObj *pShapea, const CMatrix3 &rota, const Vec3 &transla,
			GJKObj *pShapeb, const CMatrix3 &rotb, const Vec3 &translb,
			Vec3 &reta, Vec3 &retb);

class GJKSphereObj: public GJKObj
{
public:
	GJKSphereObj(float rad): mRadius(rad), mCachedPt(0, rad, 0)
	{
	}

	virtual Vec3 GetSupportPoint(const Vec3 &dir);
	virtual Vec3 GetLastPoint();
	void ResetCache();

protected:
	 float mRadius;
	 Vec3 mCachedPt;

};

class GJKCapsuleObj: public GJKObj
{
public:
	GJKCapsuleObj(float rad, float h): mRadius(rad), mHeight(h), mCachedPt(0, h + rad, 0)
	{
	}

	virtual Vec3 GetSupportPoint(const Vec3 &dir);
	virtual Vec3 GetLastPoint();
	void ResetCache();

protected:
	 float mRadius, mHeight;
	 Vec3 mCachedPt;

};




#endif