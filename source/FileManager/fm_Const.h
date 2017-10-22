#ifndef _FM_CONST_H_
#define _FM_CONST_H_

#include "all.h"

#define FM_ERR_OK			0
#define FM_ERR_NOFILE		1
#define FM_ERR_NOMEMORY		2

enum EOpenMode
{
	FM_MODE_READ,	//File will be opened for reading. File must exist.
	FM_MODE_WRITE,	//File will be created for writing, Existing file will be overwritten.
	FM_MODE_ADD		//Файл будет открыт для добавления. Если он не существует, то будет создан.
};

#endif
