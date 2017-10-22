#include "m_MathFunctions.h"

Vec3 math::NewellNormalize(const std::vector<Vec3> &perimeter)
{
	unsigned long size = perimeter.size();
	Vec3 norm = Vec3(0, 0, 0), curr, next, sum, dif;

	for(std::vector<Vec3>::const_iterator it = perimeter.begin(), end = perimeter.end(); it != end; ++it)
	{
		curr = *it;
		next = (it + 1) == perimeter.end() ? *(perimeter.begin()) : *(++it);

		sum = curr + next;
		dif = curr - next;
		norm.x += dif.y * sum.z;
		norm.y += dif.z * sum.x;
		norm.z += dif.x * sum.y;
	}

	norm.Normalize();

	return norm;
}

bool math::LooksAt(Vec3 LookPoint, Vec3 p1, Vec3 p2, Vec3 p3)
{
	return M_GREATER(Distance(LookPoint, p1, p2, p3), 0);
}

bool math::LooksAtB(Vec3 LookPoint, Vec3 p1, Vec3 p2, Vec3 p3)
{
	//-return (M_EQUAL(Dist(LookPoint, p1, p2, p3), 0) || M_GREATER(Dist(LookPoint, p1, p2, p3), 0));

	return M_GEQUAL(Distance(LookPoint, p1, p2, p3), 0);
}

FLOATTYPE math::Distance(Vec3 P, Vec3 p1, Vec3 p2, Vec3 p3)
{
	//EXPERIMENT
	Vec3 n = p2-p1;
	n = n.Cross(p3-p1);
	n.Normalize();

	return (n.Dot(P) - n.Dot(p1));
}

FLOATTYPE math::Factorial(unsigned long n)
{
	return n == 0 ? 1.0f : (FLOATTYPE)n * Factorial(n - 1);
}

FLOATTYPE math::CMN(unsigned long m, unsigned long n)
{
	if(m > n)
		return 0;

	return Factorial(n) / (Factorial(m) * Factorial(n - m));
}

FLOATTYPE math::GetLinesDistanceDOT(const Vec3 &r1, const Vec3 &r2,	const Vec3 &r3, const Vec3 &r4,
									const Vec3 &opos1, const Vec3 &opos2,
									const Vec3 &vel1, const Vec3 &vel2,
									const Vec3 &omega1, const Vec3 &omega2,
									FLOATTYPE *pDistSqrDOT, FLOATTYPE *pDist)
{
	// -- DOTs of the points of the segments --

	Vec3 r1_dt = omega1.Cross(r1 - opos1) + vel1;//Just a moving point of an object.
	Vec3 r2_dt = omega1.Cross(r2 - opos1) + vel1;

	Vec3 r3_dt = omega2.Cross(r3 - opos2) + vel2;
	Vec3 r4_dt = omega2.Cross(r4 - opos2) + vel2;

	// -- Some vectors and their DOTs --

	Vec3 d1 = r2 - r1;
	Vec3 d2 = r4 - r3;
	Vec3 dd = r1 - r3;

	Vec3 d1_dt = r2_dt - r1_dt;
	Vec3 d2_dt = r4_dt - r3_dt;
	Vec3 dd_dt = r1_dt - r3_dt;

	// -- Some constants and their DOTs --

	FLOATTYPE a = d1.Dot(d1);
	FLOATTYPE b = d1.Dot(d2);
	FLOATTYPE c = d2.Dot(d2);
	FLOATTYPE d = d1.Dot(dd);
	FLOATTYPE e = d2.Dot(dd);
	FLOATTYPE f = a * c - b * b;//a * c - b^2

	FLOATTYPE a_dt = 2.0f * d1.Dot(d1_dt);
	FLOATTYPE b_dt = d1.Dot(d2_dt) + d2.Dot(d1_dt);
	FLOATTYPE c_dt = 2.0f * d2.Dot(d2_dt);
	FLOATTYPE d_dt = d1.Dot(dd_dt) + dd.Dot(d1_dt);
	FLOATTYPE e_dt = d2.Dot(dd_dt) + dd.Dot(d2_dt);
	FLOATTYPE f_dt = a * c_dt + c * a_dt + 2.0f * b * b_dt;

	// -- The two constants and their DOTs --

	FLOATTYPE s, t;
	FLOATTYPE s_dt, t_dt;

	if(M_EQUAL_AT(f, 0.0, 0.0000001f))
	{
		s = 0.0;
		s_dt = 0.0;

		if(b > c)
		{
			t = d / b;
			t_dt = (d_dt * b - d * b_dt) / (b * b);
		}
		else
		{
			t = e / c;
			t_dt = (e_dt * c - e * c_dt) / (c * c);
		}
	}
	else
	{
		FLOATTYPE s_nom = b * e - c * d;
		FLOATTYPE t_nom = a * e - b * d;

		s = s_nom / f;
		t = t_nom / f;

		FLOATTYPE s_nom_dt = b_dt * e + b * e_dt - c_dt * d - c * d_dt;
		FLOATTYPE t_nom_dt = a_dt * e + a * e_dt - b_dt * d - b * d_dt;

		FLOATTYPE f_sqr = f * f;

		s_dt = (s_nom_dt * f - s_nom * f_dt) / f_sqr;
		t_dt = (t_nom_dt * f - t_nom * f_dt) / f_sqr;
	}

	// -- The nearest points --

	Vec3 p1 = r1 + d1 * s;
	Vec3 p2 = r2 + d2 * t;

	Vec3 p1_dt = r1_dt + d1_dt * s + d1 * s_dt;
	Vec3 p2_dt = r2_dt + d2_dt * t + d2 * t_dt;

	// -- Derivative over time of the square distance --

	FLOATTYPE dist_sqr_dt = 2.0f * ( p1.Dot(p1_dt) + p2.Dot(p2_dt) - (p1.Dot(p2_dt) + p2.Dot(p1_dt)) );

	// -- The distance --

	FLOATTYPE dist = (p2 - p1).GetLength();

	// -- Derivative over time of the distance --

	FLOATTYPE dist_dt = 0.5f * dist_sqr_dt / dist;

	// -- Other results --

	if(pDistSqrDOT)
		*pDistSqrDOT = dist_sqr_dt;

	if(pDist)
		*pDist = dist;

	return dist_dt;
}

FLOATTYPE math::GetVFDistanceDOT(const Vec3 &r, const Vec3 &v,	const Vec3 &vv,
								const Vec3 &opos1, const Vec3 &opos2,
								const Vec3 &vel1, const Vec3 &vel2,
								const Vec3 &omega1, const Vec3 &omega2)
{
	Vec3 r_dt = omega1.Cross(r - opos1) + vel1;

	Vec3 v_dt = omega2.Cross(v - opos2) + vel2;
	Vec3 vv_dt = omega2.Cross(vv - opos2) + vel2;

	Vec3 n = vv - v;
	Vec3 n_dt = vv_dt - v_dt;

	FLOATTYPE w = n.Dot(v);
	FLOATTYPE w_dt = n_dt.Dot(v) + n.Dot(v_dt);

	//FLOATTYPE d = n.Dot(r) - w;
	FLOATTYPE d_dt = n_dt.Dot(r) + n.Dot(r_dt) - w_dt;

	return d_dt;
}

Vec3 math::ChooseRandomOrthogonal(const Vec3 &arg)
{
	Vec3 ret;

	if(!M_EQUAL(arg.z, 0.0f))
	{
		ret.Set(1.0f, 1.0f, -(arg.x + arg.y) / arg.z);
	}
	else if(!M_EQUAL(arg.y, 0.0f))
	{
		ret.Set(1.0f, -(arg.x + arg.z) / arg.y, 1.0f);
	}
	else if(!M_EQUAL(arg.x, 0.0f))
	{
		ret.Set(-(arg.z + arg.y) / arg.x, 1.0f, 1.0f);
	}
	else
		ret.Set(0, 0, 0);

	return ret;
}