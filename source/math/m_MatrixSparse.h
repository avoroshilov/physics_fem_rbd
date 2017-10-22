#ifndef _M_MATRIX_SPARSE_H_
#define _M_MATRIX_SPARSE_H_

#include "all.h"

// Yale Sparse Matrix Format
class CMatrixSparse
{
public:
	CMatrixSparse(const char *, long, long);		// From Shape Matrix, 0 is empty element, other is - filled
	~CMatrixSparse(void);

	FLOATTYPE &							operator ()(unsigned long i, unsigned long j);

private:

	unsigned int						m_BlockRows;
	unsigned int						m_BlockCols;

	unsigned long						m_Rows;
	unsigned long						m_Cols;
	unsigned long						m_NNZ;		// Number of Non-Zero elements
	FLOATTYPE *							m_Matrix;

	unsigned *							m_IMatrix;
	unsigned *							m_JMatrix;
};


#endif
