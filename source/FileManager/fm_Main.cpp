#include "fm_Main.h"

#include "helpers/debug.h"
#include "helpers/log.h"

#include <winbase.h>

#include "fm_Const.h"

//using namespace fm;
using std::string;

bool fm::GetRootFolder(char *dest)
{
	assert(dest);

	char str[1024];
	//Let's get pszDestinationString.
	{
		if(GetModuleFileName(GetModuleHandle(NULL), str, 1024) == 0)
		{
			gLog.Print(LOGLEV_ERROR, "Error: in fm::GetRootFolder() - unable to get module's file name.");

			return false;
		}

		for(long c = strlen(str) - 1; c > 0; c--)
		{
			if(str[c] == '\\')
			{
				str[c] = '\0';

				break;
			}
		}

		strcat(str, "\\");

		strcpy(dest, str);
	}

	return true;
}

bool fm::GetRootFolder(string &dest)
{
	char str[1024];
	//Let's get pszDestinationString.
	{
		if(GetModuleFileName(GetModuleHandle(NULL), str, 1024) == 0)
		{
			gLog.Print(LOGLEV_ERROR, "Error: in fm::GetRootFolder() - unable to get module's file name.");

			return false;
		}

		for(long c = strlen(str) - 1; c > 0; c--)
		{
			if(str[c] == '\\')
			{
				str[c] = '\0';

				break;
			}
		}
	
		strcat(str, "\\");

		dest = str;
	}

	return true;
}

fm::CFile::CFile(void)
{
	mbIsOpened = false;
	mFile = NULL;
}

fm::CFile::~CFile(void)
{
	if(mbIsOpened)
		this->Close();
}

unsigned long fm::CFile::Open(const string &path, EOpenMode mode)
{
	assert_r(!mbIsOpened, false);

	char ModeString[10];

	switch(mode)
	{
	case FM_MODE_READ:
		strcpy(ModeString, "rb");

		break;

	case FM_MODE_WRITE:
		strcpy(ModeString, "wb");

		break;

	case FM_MODE_ADD:
		strcpy(ModeString, "ab");

		break;
	}

	mFile = fopen(path.c_str(), ModeString);

	if(!mFile)
	{
		mFile = NULL;

		return FM_ERR_NOFILE;
	}

	mbIsOpened = true;
	mMode = mode;

	return FM_ERR_OK;
}

unsigned long fm::CFile::LocalOpen(const string &path, EOpenMode mode)
{
	assert_r(!mbIsOpened, false);

	char ModeString[10];

	switch(mode)
	{
	case FM_MODE_READ:
		strcpy(ModeString, "rb");

		break;

	case FM_MODE_WRITE:
		strcpy(ModeString, "wb");

		break;

	case FM_MODE_ADD:
		strcpy(ModeString, "ab");

		break;
	}

	string root;
	GetRootFolder(root);

	string FullPath;
	FullPath.reserve(root.length() + path.length());
	FullPath.append(root);
	FullPath.append(path);

	mFile = fopen(FullPath.c_str(), ModeString);

	if(!mFile)
	{
		mFile = NULL;

		return FM_ERR_NOFILE;
	}

	mbIsOpened = true;
	mMode = mode;

	return FM_ERR_OK;
}

void fm::CFile::Close(void)
{
	if(mbIsOpened)
	{
		if(mMode == FM_MODE_WRITE || mMode == FM_MODE_ADD)
			fflush(mFile);

		fclose(mFile);

		mFile = NULL;

		mbIsOpened = false;
	}
}

unsigned long fm::CFile::GetSize(void)
{
	if(!mbIsOpened)
		return 0;

	unsigned long pos = ftell(mFile);
	fseek(mFile, 0, SEEK_END);

	unsigned long size = ftell(mFile);//Current position of the file cursor == size.

	fseek(mFile, pos, SEEK_SET);//Restore the file cursor position.

	return size;
}

unsigned long fm::CFile::GetPos(void)
{
	if(!mbIsOpened)
		return 0;

	return (unsigned long)ftell(mFile);
}

void fm::CFile::Reset(void)
{
	if(!mbIsOpened)
		return;

	fseek(mFile, 0, SEEK_SET);
}

bool fm::CFile::Seek(long bytes, int mode)
{
	if(!mbIsOpened)
		return false;

	return (fseek(mFile, bytes, mode) == 0);
}

bool fm::CFile::Read(unsigned long size, void *buffer, unsigned long *pTotallyRead)
{
	if(pTotallyRead)
		*pTotallyRead = 0;

	if(size == 0)
		return true;

	assert(mbIsOpened);
	assert(buffer);

	unsigned long ReadSize = fread(buffer, 1, size, mFile);

	if(pTotallyRead)
		*pTotallyRead = ReadSize;

	if(ReadSize < size)
		return false;

	return true;
}

bool fm::CFile::Write(unsigned long size, const void *buffer, unsigned long *pTotallyWritten)
{
	if(pTotallyWritten)
		*pTotallyWritten = 0;

	if(size == 0)
		return true;

	assert_r(mbIsOpened, false);
	assert_r(buffer, false);

	unsigned long WriteSize = fwrite(buffer, 1, size, mFile);

	if(pTotallyWritten)
		*pTotallyWritten = WriteSize;

	assert(!(WriteSize < size));//Is it not enough space or writing error?

	if(WriteSize < size)
		return false;

	return true;
}
