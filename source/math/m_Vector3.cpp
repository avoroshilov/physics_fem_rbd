#include "m_Vector3.h"

#include "helpers/log.h"

#ifndef INLINE_BUILD
	#include "m_Vector3.inl"
#endif

void Vec3::Normalize(void)
{
	FLOATTYPE magnitude = sqrt(x*x + y*y + z*z);

	if(magnitude <= M_EPSILON)
	{
		return;
	}

	x = x / magnitude;
	y = y / magnitude;
	z = z / magnitude;
}

Vec3 Vec3::GetNormalized(void) const
{
	Vec3 Result;

	FLOATTYPE magnitude = sqrt(x*x + y*y + z*z);

	if (magnitude <= M_EPSILON)
	{
		//assert(false && "Magnitude of a vector is lower than M_EPSILON!");
		gLog.Print("Magnitude of a vector is lower than M_EPSILON!");

		return Vec3(0.0, 0.0, 0.0);
	}

	Result.x = x / magnitude;
	Result.y = y / magnitude;
	Result.z = z / magnitude;

	return Result;
}

void Vec3::TangentSpace(Vec3 &Tan1, Vec3 &Tan2) const
{
	// Decide, what vector should we cross with
	if (x * x < z * z)
	{
		// Cross with vector (1 0 0)^T
		float vec_len = sqrtf(z * z + y * y);
		
		if (vec_len < M_EPSILON)
		{
			gLog.Print("Magnitude of a vector is lower than M_EPSILON [TangentSpace]!");
			return;
		}

		vec_len = 1.0f / vec_len;

		Tan1.x = 0.0f;
		Tan1.y = -vec_len * z;
		Tan1.z =  vec_len * y;
	}
	else
	{
		// Cross with vector (0 0 1)^T
		float vec_len = sqrtf(y * y + x * x);
		
		if (vec_len < M_EPSILON)
		{
			gLog.Print("Magnitude of a vector is lower than M_EPSILON [TangentSpace]!");
			return;
		}

		vec_len = 1.0f / vec_len;

		Tan1.x = -vec_len * y;
		Tan1.y =  vec_len * x;
		Tan1.z = 0.0f;
	}

	Tan2 = this->Cross(Tan1);
	Tan2.Normalize();
}
