#ifndef _SIMPLEALGS_H_
#define _SIMPLEALGS_H_

#include "math/m_Vector3.h"
#include "math/m_Matrix3.h"


float SegSeg(const Vec3 &s1a, const Vec3 &s1b, const Vec3 &s2a, const Vec3 &s2b, Vec3 &ret1, Vec3 &ret2);

float SphereSphereDistance(float rada, const Vec3 &transla,
							float radb, const Vec3 &translb,
							Vec3 &reta, Vec3 &retb);

float CapsuleSphereDistance(float rada, float ha, const CMatrix3 &rota, const Vec3 &transla,
							float radb, const Vec3 &translb,
							Vec3 &reta, Vec3 &retb);

float CapsuleCapsuleDistance(float rada, float ha, const CMatrix3 &rota, const Vec3 &transla,
							float radb, float hb, const CMatrix3 &rotb, const Vec3 &translb,
							Vec3 &reta, Vec3 &retb);

#endif