#ifndef _M_MATRIX_MN_H_
#define _M_MATRIX_MN_H_

#include "all.h"

namespace math
{
	class CMatrixMN
	{
	public:
		CMatrixMN(void);//1x1 matrix.
		CMatrixMN(unsigned long rows, unsigned long columns);

		~CMatrixMN(void);

		FLOATTYPE &							operator ()(unsigned long i, unsigned long j);
		const FLOATTYPE &						operator ()(unsigned long i, unsigned long j) const;

		//This will not copy previous contents.
		bool								SetDims(unsigned long rows, unsigned long columns);
		void								GetDims(unsigned long &rows, unsigned long &columns) const;
		unsigned long						GetRows(void) const;
		unsigned long						GetColumns(void) const;

		void								SwapRows(unsigned long i1, unsigned long i2);

	private:
		unsigned long						mRows;
		unsigned long						mColumns;
		FLOATTYPE **							mMatrix;
	};
}

#endif
