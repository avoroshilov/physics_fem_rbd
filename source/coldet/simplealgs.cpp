#include "simplealgs.h"

float SegSeg(const Vec3 &s1a, const Vec3 &s1b, const Vec3 &s2a, const Vec3 &s2b, Vec3 &ret1, Vec3 &ret2)
{
    Vec3   u = s1b - s1a;
    Vec3   v = s2b - s2a;
    Vec3   w = s1a - s2a;
	float    a = u.Dot(u);        // always >= 0
    float    b = u.Dot(v);
    float    c = v.Dot(v);        // always >= 0
    float    d = u.Dot(w);
    float    e = v.Dot(w);
    float    D = a * c - b * b;       // always >= 0
    float    sc, sN, sD = D;      // sc = sN / sD, default sD = D >= 0
    float    tc, tN, tD = D;      // tc = tN / tD, default tD = D >= 0

    // compute the line parameters of the two closest points
    if (D < M_EPSILON) 
	{ // the lines are almost parallel
        sN = 0.0;        // force using point P0 on segment S1
        sD = 1.0;        // to prevent possible division by 0.0 later
        tN = e;
        tD = c;
    }
    else
	{                // get the closest points on the infinite lines
        sN = (b*e - c*d);
        tN = (a*e - b*d);
        
		if (sN < 0.0)
		{       // sc < 0 => the s=0 edge is visible
            sN = 0.0;
            tN = e;
            tD = c;
        }
        else if (sN > sD)
		{  // sc > 1 => the s=1 edge is visible
            sN = sD;
            tN = e + b;
            tD = c;
        }
    }

    if (tN < 0.0)
	{           // tc < 0 => the t=0 edge is visible
        tN = 0.0;
        // recompute sc for this edge
        if (-d < 0.0)
		{
            sN = 0.0;
		}
        else if (-d > a)
		{
            sN = sD;
		}
        else 
		{
            sN = -d;
            sD = a;
        }
    }
    else if (tN > tD)
	{      // tc > 1 => the t=1 edge is visible
        tN = tD;
        // recompute sc for this edge
        if ((-d + b) < 0.0)
		{
            sN = 0;
		}
        else if ((-d + b) > a)
		{
            sN = sD;
		}
        else
		{
            sN = (-d + b);
            sD = a;
        }
    }

	// finally do the division to get sc and tc
    sc = (abs(sN) < M_EPSILON ? 0.0f : sN / sD);
    tc = (abs(tN) < M_EPSILON ? 0.0f : tN / tD);

    // get the difference of the two closest points
    Vec3  dP = w + (sc * u) - (tc * v);  // = S1(sc) - S2(tc)
	
	ret1 = s1a + sc * u;
	ret2 = s2a + tc * v;

	return dP.GetSquareLength();   // return the closest distance
}

float SphereSphereDistance(float rada, const Vec3 &transla,
							float radb, const Vec3 &translb,
							Vec3 &reta, Vec3 &retb)
{
	Vec3 dir = transla - translb;

	float len = dir.GetLength();
	dir.Normalize();

	reta = transla - dir * rada;
	retb = translb + dir * radb;

	float dist = len - (rada + radb);
	
	return (dist < 0 ? -1.0f : 1.0f) * dist * dist;
}

float CapsuleSphereDistance(float rada, float ha, const CMatrix3 &rota, const Vec3 &transla,
							float radb, const Vec3 &translb,
							Vec3 &reta, Vec3 &retb)
{
	Vec3 worlddir = rota * Vec3(0.0f, ha, 0.0f);
	float proj = worlddir.Dot(translb - transla);
	proj = proj < -ha * ha ? -ha * ha : proj > ha * ha ? ha * ha : proj;

	Vec3 pt = transla + worlddir * (proj / (ha * ha));
	
	Vec3 dir = pt - translb;

	float len = dir.GetLength();
	dir.Normalize();

	reta = pt - dir * rada;
	retb = translb + dir * radb;

	float dist = len - (rada + radb);
	
	return (dist < 0 ? -1.0f : 1.0f) * dist * dist;
}

float CapsuleCapsuleDistance(float rada, float ha, const CMatrix3 &rota, const Vec3 &transla,
							float radb, float hb, const CMatrix3 &rotb, const Vec3 &translb,
							Vec3 &reta, Vec3 &retb)
{
	Vec3 worlddira = rota * Vec3(0.0f, ha, 0.0f);
	Vec3 worlddirb = rotb * Vec3(0.0f, hb, 0.0f);

	Vec3 pa, pb;

	SegSeg(transla - worlddira, transla + worlddira, translb - worlddirb, translb + worlddirb, pa, pb);
	
	Vec3 d = pa - pb;

	float len = d.GetLength();
	d.Normalize();

	reta = pa - d * rada;
	retb = pb + d * radb;

	float dist = len - (rada + radb);
	
	return (dist < 0 ? -1.0f : 1.0f) * dist * dist;
}

