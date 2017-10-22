#include "all.h"
#include "gjk.h"

#include "math/m_Vector3.h"
#include "math/m_Matrix3.h"
#include "helpers/debug.h"
#include "helpers/log.h"

bool testSimplexExtraset(unsigned int testedset, unsigned int usedpoints, float dets[16][4])
{
	for(int j = 0; j < 4; ++j)
	{
		if((testedset & (1 << j)) || !(usedpoints & (1 << j)))
			continue;

		if(M_GREATER(dets[testedset | (1 << j)][j] , 0.0f))
		{
			return false;
		}
	}

	return true;
}

bool testSimplex(unsigned int testedset, unsigned int usedpoints, float dets[16][4])
{
	for(int j = 0; j < 4; ++j)
	{
		if(!(testedset & (1 << j)))
			continue;

		if(dets[testedset][j] <= 0.0f)
		{
			return false;
		}
	}

	return testSimplexExtraset(testedset, usedpoints, dets);
}

float GJK(GJKObj *pShapea, const CMatrix3 &rota, const Vec3 &transla,
			GJKObj *pShapeb, const CMatrix3 &rotb, const Vec3 &translb,
			Vec3 &reta, Vec3 &retb)
{
	const float len_tolerance = 0.1f;
	const float len_tolerance_sq = len_tolerance * len_tolerance;

	const float alg_epsilon = 0.001f;

	float simplex[4][3] = 
		{ 
			{0.0f, 0.0f, 0.0f},
			{0.0f, 0.0f, 0.0f},
			{0.0f, 0.0f, 0.0f},
			{0.0f, 0.0f, 0.0f}
		}; //at most four points, xyz each
	
	Vec3 pointsa[4], pointsb[4];

	float dots[4][4] = 
		{
			{0.0f, 0.0f, 0.0f, 0.0f},
			{0.0f, 0.0f, 0.0f, 0.0f},
			{0.0f, 0.0f, 0.0f, 0.0f},
			{0.0f, 0.0f, 0.0f, 0.0f}
		}; //dot products of the points on each other
	
	float dets[16][4] = 
		{
			{0.0f, 0.0f, 0.0f, 0.0f}, //not used
			{1.0f, 0.0f, 0.0f, 0.0f}, //only 1 used
			{0.0f, 1.0f, 0.0f, 0.0f}, //only 2 used
			{0.0f, 0.0f, 0.0f, 0.0f}, //only 1 and 2 are used

			{0.0f, 0.0f, 1.0f, 0.0f}, //only 3 is used
			{0.0f, 0.0f, 0.0f, 0.0f}, //only 1 and 3 used
			{0.0f, 0.0f, 0.0f, 0.0f}, //only 2 and 3 used
			{0.0f, 0.0f, 0.0f, 0.0f}, //only 1 2 3 used

			{0.0f, 0.0f, 0.0f, 1.0f}, //only 4 used
			{0.0f, 0.0f, 0.0f, 0.0f}, //1 and 4
			{0.0f, 0.0f, 0.0f, 0.0f}, //2 and 4
			{0.0f, 0.0f, 0.0f, 0.0f}, //1, 2 and 4

			{0.0f, 0.0f, 0.0f, 0.0f}, //3 and 4
			{0.0f, 0.0f, 0.0f, 0.0f}, //1, 3 and 4
			{0.0f, 0.0f, 0.0f, 0.0f}, //1, 2 and 4
			{0.0f, 0.0f, 0.0f, 0.0f}  //1, 2, 3 and 4
		}; //delta-i determinant; first index is bit-array simplex descriptor

	//if we get result at once, use the first point
	float lambdas[4] = {1.0f, 0.0f, 0.0f, 0.0f};

	unsigned int usedpoints = 0, oldusedpoints = 0;

	//arbitary point in A-B space
	pointsa[0] = rota * pShapea->GetLastPoint() + transla;
	pointsb[0] = (rotb * pShapeb->GetLastPoint() + translb);
	Vec3 v = pointsa[0] - pointsb[0]; 
	pointsa[1] = pointsb[1] = pointsa[2] = pointsb[2] = pointsa[3] = pointsb[3] = Vec3(0.0f, 0.0f, 0.0f);
	bool bContinue = true;
	unsigned int const maxiter = 300;
	unsigned int i;
	float mu = 0.0f;
	Vec3 support, supporta, supportb;

	gLog2.LogDebug("\n\nNew call\n");

	for(i = 0; bContinue && i < maxiter; ++i)
	{
		float v_sq_len = v.GetSquareLength();
	
		gLog2.LogDebug("%f %f %f : %f\n", v.x, v.y, v.z, v_sq_len);
	
#ifdef DEBUG_BUILD
		
		Vec3 simpoints[4];
		int simpoint = 0;

		for(int j = 0; j < 4; ++j)
		{
			if(!(usedpoints & (1 << j)))
			{
				continue;
			}
			
			gLog2.LogDebug("Simplex point %d: %f %f %f\n", j, simplex[j][0], simplex[j][1], simplex[j][2]); 
			simpoints[simpoint++] = Vec3(simplex[j][0], simplex[j][1], simplex[j][2]);  
		}

		float simplexsize = simpoint == 0 ? 0 : simpoint == 1 ? 1.0f : simpoint == 2 ? (simpoints[1] - simpoints[0]).GetLength() :
			simpoint == 3 ? ((simpoints[2] - simpoints[0]).Cross(simpoints[1] - simpoints[0])).GetSquareLength() : 666.0f;
		
		gLog2.LogDebug("Simplex size %f ; %d simplex points\n", simplexsize, simpoint);

#endif

		if(v_sq_len < len_tolerance_sq)
		{
			bContinue = false;
			break;
		}
		
		//support point of A-B in the direction -v
		Vec3 vnorm = v.GetNormalized();
		supporta = transla + rota * pShapea->GetSupportPoint(-vnorm * rota);
		supportb = translb + rotb * pShapeb->GetSupportPoint(vnorm * rotb);
		support = supporta - supportb;
		
		float sigma = v.Dot(support);
		mu = mu > sigma ? mu : sigma;

		if(v_sq_len - mu <= alg_epsilon * v_sq_len)
		{
			bContinue = false;
			break;
		}

		//test if we arrive at the same point - degenerate check

		unsigned int newslot = 0xFF;

		for(int j = 0; j < 4; ++j)
		{
			if(!(usedpoints & (1 << j)))
			{
				newslot = j;
				break;
			}
		}
	
		for(int j = 0; j < 4; ++j)
		{
			if(!(oldusedpoints & (1 << j)))
			{
				continue;
			}
			
			if(M_EQUAL(simplex[j][0], support.x) && M_EQUAL(simplex[j][1], support.y) && M_EQUAL(simplex[j][2], support.z))
			{
				bContinue = false;
				break;
			}
		}

		if(!bContinue)
			break;

		assert_r(newslot != 0xFF, -1.0f);
		
		simplex[newslot][0] = support.x;
		simplex[newslot][1] = support.y;
		simplex[newslot][2] = support.z;
		
		pointsa[newslot] = supporta;
		pointsb[newslot] = supportb;

		usedpoints |= (1 << newslot);
		oldusedpoints = usedpoints;
		
		//recalculate dots
		for(unsigned int j = 0; j < newslot; ++j)
		{
			if(!(usedpoints & (1 << j)))
				continue;

			dots[j][newslot] = dots[newslot][j] = simplex[j][0] * support.x + simplex[j][1] * support.y
																			+ simplex[j][2] * support.z;
		}
		
		dots[newslot][newslot] = support.Dot(support);

		for(int j = newslot + 1; j < 4; ++j)
		{
			if(!(usedpoints & (1 << j)))
				continue;

			dots[j][newslot] = dots[newslot][j] = simplex[j][0] * support.x + simplex[j][1] * support.y
																			+ simplex[j][2] * support.z;
		}
		
		//recalculate dets
		unsigned int newslotbit = 1 << newslot;

		//pairs
		for(int j = 0; j < 4; ++j)
		{
			if(j == newslot)
				continue;

			if(!(usedpoints & (1 << j)))
				continue;

			unsigned int bits = (1 << j) | newslotbit;
			
			dets[bits][newslot] = dots[j][j] - dots[j][newslot];
			dets[bits][j] = dots[newslot][newslot] - dots[newslot][j];
		}

		//tripples
		unsigned int shifts[3][2] = {{1, 2}, {1, 3}, {2, 3}};

		for(int j = 0; j < 3; ++j)
		{
			unsigned int ind1 = newslot + shifts[j][0];
			ind1 -= ind1 > 3 ? 4 : 0;
			unsigned int ind2 = newslot + shifts[j][1];
			ind2 -= ind2 > 3 ? 4 : 0;
					
			unsigned int bits12 = (1 << ind1) | (1 << ind2);
			
			if((usedpoints & bits12) != bits12)
				continue;

			unsigned int bits02 = newslotbit | (1 << ind2);
			unsigned int bits01 = (1 << ind1) | newslotbit;
			unsigned int bits = (1 << ind1) | (1 << ind2) | newslotbit;
			
			dets[bits][newslot] = dets[bits12][ind1] * (dots[ind1][ind1] - dots[ind1][newslot])
								+ dets[bits12][ind2] * (dots[ind2][ind1] - dots[ind2][newslot]);
			dets[bits][ind1] = dets[bits02][newslot] * (dots[newslot][newslot] - dots[newslot][ind1])
								+ dets[bits02][ind2] * (dots[ind2][newslot] - dots[ind2][ind1]);
			dets[bits][ind2] = dets[bits01][newslot] * (dots[newslot][newslot] - dots[newslot][ind2])
								+ dets[bits01][ind1] * (dots[ind1][newslot] - dots[ind1][ind2]);

		}

		//quadripple
		if(usedpoints == 0xF)
		{
			for(int j = 0; j < 4; ++j)
			{
				unsigned int bitsnoj = ~(1 << j) & 0xF;
		
				dets[0xf][j] = 0;
				
				unsigned int fixedind = j == 0 ? 1 : 0;

				for(int k = 0; k < 4; ++k)
				{
					if(k == j)
						continue;

					dets[0xf][j] += dets[bitsnoj][k] * (dots[k][fixedind] - dots[k][j]);
				}
			}
		}

		//now we must find a new simplex
		
		bool bFound = false;
		unsigned int testedset = 0;

		unsigned int maxdim = ((usedpoints & 0x8) ? 1 : 0) + ((usedpoints & 0x4) ? 1 : 0) + ((usedpoints & 0x2) ? 1 : 0)
							+ ((usedpoints & 0x1) ? 1 : 0);

		testedset = (1 << newslot);

		//for 1-d, first condition must be trivial
		bFound = testSimplexExtraset(testedset, usedpoints, dets);
		
		//test 2-d
		if(!bFound && maxdim > 1)
		{
			for(int a = 0; a < 4; ++a)
			{
				if(a == newslot || !(usedpoints & (1 << a)))
					continue;

				testedset = (1 << a) | (1 << newslot);
				
				if(testSimplex(testedset, usedpoints, dets))
				{
					bFound = true;
					break;
				}
			}
		}

		//test 3-d
		if(!bFound && maxdim > 2)
		{
			for(int a = 0; a < 4; ++a)
			{
				if(a == newslot || !(usedpoints & (1 << a)))
					continue;

				for(int b = a + 1; b < 4; ++b)
				{
					if(b == newslot || !(usedpoints & (1 << b)))
						continue;
					
					testedset = (1 << a) | (1 << b) | (1 << newslot);
					
					if(testSimplex(testedset, usedpoints, dets))
					{
						bFound = true;
						break;
					}
				}

				if(bFound)
					break;
			}
		}

		//test 4-d
		if(!bFound && maxdim > 3)
		{
			testedset = 0xF;
			
			if(testSimplex(testedset, usedpoints, dets))
			{
				bFound = true;
				break;
			}
		}
		
		if(!bFound)
		{
			//degenerate case
			bContinue = false;
			break;
		}

		usedpoints = testedset;

		float bigdet = 0.0;

		for(int j = 0; j < 4; ++j)
		{
			if(testedset & (1 << j))
			{
				bigdet += dets[testedset][j];
			}
		}
		
		for(int j = 0; j < 4; ++j)
		{
			if(testedset & (1 << j))
			{
				lambdas[j] = dets[testedset][j] / bigdet;
			}
			else
			{
				lambdas[j] = 0.0f;
			}
		}
	
		v = lambdas[0] * Vec3(simplex[0][0], simplex[0][1], simplex[0][2]) + 
			lambdas[1] * Vec3(simplex[1][0], simplex[1][1], simplex[1][2]) +
			lambdas[2] * Vec3(simplex[2][0], simplex[2][1], simplex[2][2]) +
			lambdas[3] * Vec3(simplex[3][0], simplex[3][1], simplex[3][2]);
		
	}

	if(i >= maxiter)
		_asm{int 3};

	assert(i < maxiter); //this crap must never ever cycle forever
	
	if(M_EQUAL(lambdas[0], 0.0f) && M_EQUAL(lambdas[1], 0.0f) && M_EQUAL(lambdas[2], 0.0f) && M_EQUAL(lambdas[3], 0.0f))
		return -1.0f;

	reta = lambdas[0] * pointsa[0] + lambdas[1] * pointsa[1] + lambdas[2] * pointsa[2] + lambdas[3] * pointsa[3];
	retb = lambdas[0] * pointsb[0] + lambdas[1] * pointsb[1] + lambdas[2] * pointsb[2] + lambdas[3] * pointsb[3];

	return (retb - reta).GetSquareLength();
}

Vec3 GJKSphereObj::GetSupportPoint(const Vec3 &dir)
{
	return mCachedPt = mRadius * dir;
}

Vec3 GJKSphereObj::GetLastPoint()
{
	return mCachedPt;
}

void GJKSphereObj::ResetCache()
{
	mCachedPt.Set(0.0f, mRadius, 0.0f);
}

Vec3 GJKCapsuleObj::GetSupportPoint(const Vec3 &dir)
{
	return mCachedPt = mRadius * dir + Vec3(0.0f, dir.y > 0.0f ? mHeight : -mHeight, 0.0f);
}

Vec3 GJKCapsuleObj::GetLastPoint()
{
	return mCachedPt;
}

void GJKCapsuleObj::ResetCache()
{
	mCachedPt.Set(0.0f, mRadius + mHeight , 0.0f);
}
