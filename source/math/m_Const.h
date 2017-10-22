#ifndef _M_CONST_H_
#define _M_CONST_H_

#include "all.h"

#ifdef M_PI
	#undef M_PI
	#define M_PI (3.14159265358979323846264338327950288419716939931148196659300057)
#else
	#define M_PI (3.14159265358979323846264338327950288419716939931148196659300057f)
#endif

#define M_EPSILON (0.0001f)

#define M_ABS(x) (((x) > 0) ? (x) : (-(x)))

#define M_EQUAL(x, y)			(M_ABS((x) - (y)) <= M_EPSILON)		//x == y ? true : false;
#define M_LESS(x, y)			((x) <  ((y) - M_EPSILON))			//x < y  ? true : false;
#define M_GREATER(x, y)			((x) >  ((y) + M_EPSILON))			//x > y  ? true : false;
#define M_LEQUAL(x, y)			((x) <= ((y) + M_EPSILON))			//x <= y ? true : false;
#define M_GEQUAL(x, y)			((x) >= ((y) - M_EPSILON))			//x >= y ? true : false;

#define M_EQUAL_AT(x, y, e)		(M_ABS((x) - (y)) <= e)				//x == y ? true : false;
#define M_LESS_AT(x, y, e)		((x) <  ((y) - e))					//x < y  ? true : false;
#define M_GREATER_AT(x, y, e)	((x) >  ((y) + e))					//x > y  ? true : false;
#define M_LEQUAL_AT(x, y, e)	((x) <= ((y) + e))					//x <= y ? true : false;
#define M_GEQUAL_AT(x, y, e)	((x) >= ((y) - e))					//x >= y ? true : false;

#endif
