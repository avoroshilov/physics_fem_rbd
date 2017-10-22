#include "log.h"

#include "debug.h"

using namespace fm;

#define NO_LOGGING

/////////////////////////////////////////////////////////////////////////////////////
//class CLogChannel
/////////////////////////////////////////////////////////////////////////////////////

CLogChannel::CLogChannel(const char *LogFileName)
{
	mbInitialized = Initialize(LogFileName);
	mbEnabled = true;
}

CLogChannel::~CLogChannel(void)
{
	if(mbInitialized)
		Release();
}

void CLogChannel::SetEnabled(bool enabled)
{
	mbEnabled = enabled;
}

bool CLogChannel::GetEnabled() const
{
	return mbEnabled;
}

bool CLogChannel::Initialize(const char *LogFileName)
{
	if(mbInitialized)
		Release();

	assert_r(LogFileName, false);

	strcpy(mLogFileName, LogFileName);

	if(mLogFile.LocalOpen(mLogFileName, FM_MODE_WRITE) != FM_ERR_OK)
	{
		assert_ex(false, "CLogChannel::Initialize(const char *LogFileName) -> if(mLogFile.Open(mLogFileName, WRITE) != FM_ERR_OK) failed", false);

		return false;
	}

	mLogFile.Close();

	mDepth = 0;
	mLogLevel = LOGLEV_NOLOG;// LOGLEV_MESSAGE;
	mbInitialized = true;

	return true;
}

void CLogChannel::Release(void)
{
	assert_vr(mbInitialized);
}

void CLogChannel::SetLevel(ELogLevel loglevel)
{
	mLogLevel = loglevel;
}

ELogLevel CLogChannel::GetLevel(void)
{
	return mLogLevel;
}

void CLogChannel::IncDepth(void)
{
#ifdef NO_LOGGING
	return;
#endif

	assert_ex(mbInitialized, "Use of uninitialized CLogChannel class instance", false);

	Put("{\r\n");

	mDepth++;
}

void CLogChannel::IncDepth(const char *DepthName)
{
#ifdef NO_LOGGING
	return;
#endif

	assert_ex(mbInitialized, "Use of uninitialized CLogChannel class instance", false);
	assert_vr(DepthName);

	char str[4096];
	
	if(strlen(DepthName) + 11 > 4095)
	{
		Put("{//NAME OUTPUT FAILED.\r\n");
	}
	else
	{
		sprintf(str, "{//%s.\r\n", DepthName);
		Put(str);
	}

	mDepth++;
}

void CLogChannel::DecDepth(void)
{
#ifdef NO_LOGGING
	return;
#endif

	assert_ex(mbInitialized, "Use of uninitialized CLogChannel class instance", false);

	if(mDepth > 0)
		mDepth--;

	Put("}\r\n");
}

void CLogChannel::Put(const char *string)
{
#ifdef NO_LOGGING
	return;
#endif
	if (!mbEnabled)
		return;

	assert_ex(mbInitialized, "Use of uninitialized CLogChannel class instance", false);
	assert_vr(string);

	if(mLogFile.LocalOpen(mLogFileName, FM_MODE_ADD) == FM_ERR_OK)
	{
		for(unsigned long h = 0; h < mDepth; h++)
			mLogFile.Write(1, "\t");

		mLogFile.Write(strlen(string), string);
		mLogFile.Write(2, "\r\n");
		mLogFile.Close();
	}
}

void CLogChannel::Print(const char *message, ...)
{
#ifdef NO_LOGGING
	return;
#endif
	if (!mbEnabled)
		return;

	assert_ex(mbInitialized, "Use of uninitialized CLogChannel class instance", false);

	char mstr[10240];

	va_list argptr;
	va_start(argptr, message);
	vsprintf(mstr, message, argptr);
	va_end(argptr);

	Put(mstr);
}

void CLogChannel::Print(ELogLevel level, const char *message)
{
#ifdef NO_LOGGING
	return;
#endif
	if (!mbEnabled)
		return;

	assert_ex(mbInitialized, "Use of uninitialized CLogChannel class instance", false);

	if(level > mLogLevel || level == LOGLEV_NOLOG)
		return;
	
	char string[10240];
	static const char *logLevels[] = {
		"Nolog",
		"Crit",
		"Err",
		"Warn",
		"Msg",
		"Vrbs",
		"Dbg"
	};

	sprintf(string, "<%s> %s", logLevels[level], message);

	Put(string);
}

void CLogChannel::LogCritical(const char *message, ...)
{
#ifdef NO_LOGGING
	return;
#endif
	if (!mbEnabled)
		return;

	va_list argptr;
	va_start(argptr, message);
	char mstr[10240];
	vsprintf(mstr, message, argptr);
	va_end(argptr);
	Print(LOGLEV_CRITICAL, mstr);
}

void CLogChannel::LogError(const char *message, ...)
{
#ifdef NO_LOGGING
	return;
#endif
	if (!mbEnabled)
		return;

	va_list argptr;
	va_start(argptr, message);
	char mstr[10240];
	vsprintf(mstr, message, argptr);
	va_end(argptr);
	Print(LOGLEV_ERROR, mstr);
}

void CLogChannel::LogWarning(const char *message, ...)
{
#ifdef NO_LOGGING
	return;
#endif
	if (!mbEnabled)
		return;

	va_list argptr;
	va_start(argptr, message);
	char mstr[10240];
	vsprintf(mstr, message, argptr);
	va_end(argptr);
	Print(LOGLEV_WARNING, mstr);
}

void CLogChannel::LogMessage(const char *message, ...)
{
#ifdef NO_LOGGING
	return;
#endif
	if (!mbEnabled)
		return;

	va_list argptr;
	va_start(argptr, message);
	char mstr[10240];
	vsprintf(mstr, message, argptr);
	va_end(argptr);
	Print(LOGLEV_MESSAGE, mstr);
}

void CLogChannel::LogVerbose(const char *message, ...)
{
#ifdef NO_LOGGING
	return;
#endif
	if (!mbEnabled)
		return;

	va_list argptr;
	va_start(argptr, message);
	char mstr[10240];
	vsprintf(mstr, message, argptr);
	va_end(argptr);
	Print(LOGLEV_VERBOSE, mstr);
}

void CLogChannel::LogDebug(const char *message, ...)
{
#ifdef NO_LOGGING
	return;
#endif
	if (!mbEnabled)
		return;

	va_list argptr;
	va_start(argptr, message);
	char mstr[10240];
	vsprintf(mstr, message, argptr);
	va_end(argptr);
	Print(LOGLEV_DEBUG, mstr);
}


/////////////////////////////////////////////////////////////////////////////////////
//class CLogSection
/////////////////////////////////////////////////////////////////////////////////////

CLogSection::CLogSection(CLogChannel &channel, ELogLevel sectionlevel): mMyChannel(channel)
{
	mOldLogLevel = channel.GetLevel();
	channel.SetLevel(sectionlevel);
	mMyDepth = sNestedDepth++;
}

CLogSection::~CLogSection(void)
{
	assert(--sNestedDepth == mMyDepth && "Debug log sections overlap!");
	mMyChannel.SetLevel(mOldLogLevel);
}

unsigned int CLogSection::sNestedDepth = 0;
