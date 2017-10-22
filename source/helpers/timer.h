#ifndef _HELPERS_TIMER_H_
#define _HELPERS_TIMER_H_

#include <Windows.h>

class Timer
{
public:

	Timer()
	{
		TimerFrequency = GetFrequency();
	}

	void Start();
	double Time();

private:

	__int64 Count();

	static inline __int64 CalculateFrequency();
	static __int64 GetFrequency();

	__int64 StartTime;
	double ElapsedTime;
	__int64 TimerFrequency;
};

__int64 Timer::CalculateFrequency()
{
	LARGE_INTEGER tFreq;
	QueryPerformanceFrequency(&tFreq);
	return tFreq.QuadPart;
}

#endif