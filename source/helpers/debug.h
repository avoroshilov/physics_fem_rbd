#ifndef _DEBUG_H_
#define _DEBUG_H_

#include "assert.h"
#include "all.h"

void message(const char *str, ...);

namespace com
{
	__int64 GetCRC(void *data, unsigned long size);
}

#ifdef DEBUG_BUILD
	#include <windows.h>
	#include <stdio.h>

	//assert_t macro always terminates the process.
	#define assert_t(a)\
		if(!(a))\
		{\
			char filesrc[256];\
			strcpy(filesrc, __FILE__);\
			char file[256];\
			for(int i = strlen(filesrc) - 1; i >= 0; i--)\
			{\
				if(filesrc[i] == '\\')\
				{\
					strcpy(file, filesrc + i + 1);\
					break;\
				}\
			}\
			char s[1024];\
			sprintf(s, "\"%s\" condition falied.\r\nFile: %s\r\nLine: %ld", #a, file, __LINE__);\
			MessageBox(0, s, "Debug error", MB_ICONHAND | MB_SYSTEMMODAL);\
			TerminateProcess(GetCurrentProcess(), 0);\
		}

	//assert_ex allows additional message and terminates the process upon your choice.
	#define assert_ex(condition, message, terminate)\
		if(!(condition))\
		{\
			char s[1024];\
			sprintf(s, "\"%s\" condition falied.\r\nFile: %s\r\nLine: %ld.\r\n\r\n\
Additional message: %s", #condition, __FILE__, __LINE__, message);\
			MessageBox(0, s, "Debug error", MB_ICONHAND | MB_SYSTEMMODAL);\
			if(terminate)\
				TerminateProcess(GetCurrentProcess(), 0);\
		}

	//assert_t returns the retval if condition is wrong.
	#define assert_r(a, retval)\
		if(!(a))\
		{\
			char filesrc[256];\
			strcpy(filesrc, __FILE__);\
			char file[256];\
			for(int i = strlen(filesrc) - 1; i >= 0; i--)\
			{\
				if(filesrc[i] == '\\')\
				{\
					strcpy(file, filesrc + i + 1);\
					break;\
				}\
			}\
			char s[1024];\
			sprintf(s, "\"%s\" condition falied.\r\nFile: %s\r\nLine: %ld", #a, file, __LINE__);\
			MessageBox(0, s, "Debug error", MB_ICONHAND | MB_SYSTEMMODAL);\
			return (retval);\
		}

	//assert_vr returns no value by operator return if condition is wrong.
	#define assert_vr(a)\
		if(!(a))\
		{\
			char filesrc[256];\
			strcpy(filesrc, __FILE__);\
			char file[256];\
			for(int i = strlen(filesrc) - 1; i >= 0; i--)\
			{\
				if(filesrc[i] == '\\')\
				{\
					strcpy(file, filesrc + i + 1);\
					break;\
				}\
			}\
			char s[1024];\
			sprintf(s, "\"%s\" condition falied.\r\nFile: %s\r\nLine: %ld", #a, file, __LINE__);\
			MessageBox(0, s, "Debug error", MB_ICONHAND | MB_SYSTEMMODAL);\
			return;\
		}
#else
	#define assert_t(a)
	#define assert_ex(condition, message, terminate)
	#define assert_r(a, retval) if(!(a)) return (retval);
	#define assert_vr(a) if(!(a)) return;
#endif

#endif
