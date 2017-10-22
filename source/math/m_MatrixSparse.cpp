#include "m_MatrixSparse.h"

CMatrixSparse::CMatrixSparse(const char *ShapeMatrix, long Rows, long Cols)
{
	unsigned i, j;
	m_NNZ = 0;

	m_Rows = Rows;
	m_Cols = Cols;

	long RowShift = 0;
	for (i = 0; i < m_Rows; ++i)
	{
		for (j = 0; j < m_Cols; ++j)
			if (ShapeMatrix[RowShift + j])
			{
				++m_NNZ;
			}

			RowShift += m_Cols;
	}

	m_Matrix	= new FLOATTYPE[m_NNZ * (m_BlockRows * m_BlockCols)];
	m_IMatrix	= new unsigned[m_Rows + 1];
	m_JMatrix	= new unsigned[m_NNZ];
}


CMatrixSparse::~CMatrixSparse()
{
	if (m_Matrix)	delete [] m_Matrix;
	if (m_IMatrix)	delete [] m_IMatrix;
	if (m_JMatrix)	delete [] m_JMatrix;
}
