#ifndef _FM_IFILE_H_
#define _FM_IFILE_H_

#include "fm_Const.h"

#include <stdio.h>
#include <wtypes.h>
#include <string>


namespace fm
{
	//File access interface

	class IFile
	{
	public:
		virtual ~IFile(void)
		{
		};

		//Opens file by global path.
		virtual unsigned long	Open(const std::string &path, EOpenMode mode) = 0;

		//Opens file by relative to the root folder path.
		virtual unsigned long	LocalOpen(const std::string &path, EOpenMode mode) = 0;

		virtual void			Close(void) = 0;

		virtual unsigned long	GetSize(void) = 0;
		virtual unsigned long	GetPos(void) = 0;
		virtual void			Reset(void) = 0;
		virtual bool			Seek(long bytes, int mode) = 0;//mode == {SEEK_*}//MAKE VERIFIER FOR SEEKED BYTES!

		virtual bool			Read(unsigned long size, void *buffer, unsigned long *pTotallyRead = NULL) = 0;
		virtual bool			Write(unsigned long size, const void *buffer, unsigned long *pTotallyWritten = NULL) = 0;

	};
}

#endif
