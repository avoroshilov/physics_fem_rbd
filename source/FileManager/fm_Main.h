#ifndef _FM_MAIN_H_
#define _FM_MAIN_H_

/*****************************************************************************\
  Путь, указываемый при открытии файла - относителен. Он будет добавлен к
  рабочей папке программы.

  Файлы открываются сначала из папок, затем из паков. Файлы сохраняются
  только в папки.
\*****************************************************************************/

#include "all.h"

#include "fm_IFile.h"

#include "fm_Const.h"

#include <stdio.h>
#include <wtypes.h>
#include <string>

namespace fm
{
	bool GetRootFolder(char *dest);
	bool GetRootFolder(std::string &dest);

	class CFile: public IFile
	{
	public:
		CFile(void);
		~CFile(void);

		//Opens file by global path.
		unsigned long			Open(const std::string &path, EOpenMode mode);

		//Opens file by relative to the root folder path.
		unsigned long			LocalOpen(const std::string &path, EOpenMode mode);

		void					Close(void);

		unsigned long			GetSize(void);
		unsigned long			GetPos(void);
		void					Reset(void);
		bool					Seek(long bytes, int mode);//mode == {SEEK_*}//MAKE VERIFIER FOR SEEKED BYTES!

		bool					Read(unsigned long size, void *buffer, unsigned long *pTotallyRead = NULL);
		bool					Write(unsigned long size, const void *buffer, unsigned long *pTotallyWritten = NULL);

	private:
		bool					mbIsOpened;
		EOpenMode				mMode;
		FILE					*mFile;
	};
}

#endif
