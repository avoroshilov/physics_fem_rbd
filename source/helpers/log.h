#ifndef _LOG_H_
#define _LOG_H_

#include "all.h"

#include "FileManager/fm_Main.h"

enum ELogLevel
{
	LOGLEV_NOLOG,		//no message shall have this level
	LOGLEV_CRITICAL,	//critical error - the program is expected to be terminated
	LOGLEV_ERROR,		//non-critical error - the program is expected to behave incorrectly
	LOGLEV_WARNING,		//undesired unexpected behavior - program is expected to recover and go on
	LOGLEV_MESSAGE,		//message to be shown always in course of normal execution of the program
	LOGLEV_VERBOSE,		//extra information to be shown in corse of normal execution on demand
	LOGLEV_DEBUG		//information intersting for debugging only. May rwsult in considerable slowdown
};

class CLogChannel
{
public:

	CLogChannel(const char *LogFileName);
	~CLogChannel(void);

	void									SetEnabled(bool enabled);
	bool									GetEnabled() const;
	bool									Initialize(const char *LogFileName);
	void									Release(void);

	void									SetLevel(ELogLevel loglevel);
	ELogLevel								GetLevel(void);

	void									IncDepth(void);
	void									IncDepth(const char *DepthName);
	void									DecDepth(void);

	void									Put(const char *string);
	void									Print(const char *message, ...);
	void									Print(ELogLevel level, const char *message);

	void									LogCritical(const char *message, ...);
	void									LogError(const char *message, ...);
	void									LogWarning(const char *message, ...);
	void									LogMessage(const char *message, ...);
	void									LogVerbose(const char *message, ...);
	void									LogDebug(const char *message, ...);

protected:

	bool									mbEnabled;
	bool									mbInitialized;

	char									mLogFileName[256];
	fm::CFile								mLogFile;
	ELogLevel								mLogLevel;
	unsigned long							mDepth;
	
};

class CLogSection
{
public:

	explicit CLogSection(CLogChannel &channel, ELogLevel sectionlevel);
	~CLogSection(void);

private:
	ELogLevel						mOldLogLevel;
	unsigned int					mMyDepth;
	CLogChannel						&mMyChannel;

	static unsigned int				sNestedDepth;
};

#include "logchannels.inl"

#endif
