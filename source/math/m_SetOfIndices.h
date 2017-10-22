#ifndef _M_SET_OF_INDICES_H_
#define _M_SET_OF_INDICES_H_

#include "all.h"

namespace math
{
	class CSetOfIndices
	{
	public:
		explicit CSetOfIndices(unsigned long MaxSize);
		~CSetOfIndices(void);

		unsigned long &						operator [](unsigned long i);
		const unsigned long &				operator [](unsigned long i) const;

		unsigned long						GetSize(void) const {return mSize;}

		void								Add(unsigned long index);
		void								Delete(unsigned long index);
		void								Empty(void);

		bool								Contains(unsigned long index);

	private:
		unsigned long						mSize;
		unsigned long						mMaxSize;
		unsigned long *						mIndices;
	};
}

#endif
