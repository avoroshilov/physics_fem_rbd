#ifndef _ALL_H_
#define _ALL_H_

//#define _CRTDBG_MAP_ALLOC 1
#include <stdlib.h>
//#include <crtdbg.h>

#ifdef DEBUG_BUILD
	#error DEBUG_BUILD is defined!
#endif

//Comment this line, if you want to get exe without debug parts of code.
#define DEBUG_BUILD

#ifdef INLINE_BUILD
	#error INLINE_BUILD is defined!
#endif

#ifdef INLINE
	#error INLINE is defined!
#endif

//Comment this line, if you want to get exe without inlining functions.
#define INLINE_BUILD

#ifdef INLINE_BUILD
	#define INLINE __inline
#else
	#define INLINE
#endif

#define WIN32_LEAN_AND_MEAN

#pragma warning (disable: 4996) //strcpy declared deprecated
#pragma warning (disable: 4267) //strcpy declared deprecated

//#pragma warning (disable: 4244) //conversions, possible loss of data
//#pragma warning (disable: 4305) //truncations, possible loss of data

//#pragma warning (disable: 1572) // floating-point ==

//#define for if(true)for //To make for's variable its local variable.

#define DIRECTINPUT_VERSION 0x0800

#define FLOATTYPE float
#define FLOATSUFFIX f

/*inline void* __cdecl operator new(unsigned int s)
{ 
	void* p =  malloc(s);
	if(s == 52)
		char b = 'b';
	return p; 
};*/

#include "helpers/globalallocator.h"

#endif

