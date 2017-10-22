#include "debug.h"

void message(const char *str, ...)
{
	#ifdef DEBUG_BUILD

	char info[10240];
	va_list argptr;
	va_start(argptr, str);
	vsprintf(info, str, argptr);
	va_end(argptr);

	MessageBox(0, info, "Debug message", MB_ICONINFORMATION | MB_SYSTEMMODAL);

	#endif
}

__int64 com::GetCRC(void *data, unsigned long size)
{
	assert_t(data);

	unsigned char *p = (unsigned char *)data;
	__int64 res = 0;

	for(unsigned long i = 0; i < size; i++)
		res += p[i];

	return res;
}
