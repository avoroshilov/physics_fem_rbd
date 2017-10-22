#include "timer.h"

void Timer::Start(void)
{
	LARGE_INTEGER s;
	QueryPerformanceCounter(&s);
	StartTime = s.QuadPart;
}

double Timer::Time(void)
{
	ElapsedTime = ((Count() - StartTime) * 1000 / (double)TimerFrequency);
	return ElapsedTime;
}

//////////////////////////////////////////////////////////////////////////
// Private functions
//////////////////////////////////////////////////////////////////////////

__int64 Timer::Count(void)
{
	LARGE_INTEGER s;
	QueryPerformanceCounter(&s);
	return s.QuadPart;
}

__int64 Timer::GetFrequency()
{
	static __int64 Freq = CalculateFrequency();
	return Freq;
}

