#include "m_MatrixMN.h"

#include "helpers/debug.h"

using namespace math;

CMatrixMN::CMatrixMN(void)
{
	mMatrix = new FLOATTYPE *[1];
	mMatrix[0] = new FLOATTYPE[1];

	mRows = 1;
	mColumns = 1;
}

CMatrixMN::CMatrixMN(unsigned long rows, unsigned long columns)
{
	assert(rows > 0 && columns > 0);

	if(rows == 0)
		rows = 1;

	if(columns == 0)
		columns = 1;

	mMatrix = new FLOATTYPE *[rows];

	assert_t(mMatrix && "Error: not enough memory!");

	for(unsigned long i = 0; i < rows; i++)
	{
		mMatrix[i] = new FLOATTYPE[columns];
		assert_t(mMatrix[i] && "Error: not enough memory!");
	}

	mRows = rows;
	mColumns = columns;
}

CMatrixMN::~CMatrixMN(void)
{
	if(mMatrix)
	{
		for(unsigned long i = 0; i < mRows; i++)
			delete [] mMatrix[i];

		delete [] mMatrix;

		mMatrix = NULL;
	}
}

FLOATTYPE &CMatrixMN::operator ()(unsigned long i, unsigned long j)
{
	assert_t(i < mRows && j < mColumns && "Index range check error.");

	return mMatrix[i][j];
}

const FLOATTYPE &CMatrixMN::operator ()(unsigned long i, unsigned long j) const
{
	assert_t(i < mRows && j < mColumns && "Index range check error.");

	return mMatrix[i][j];
}

bool CMatrixMN::SetDims(unsigned long rows, unsigned long columns)
{
	assert(rows > 0 && columns > 0);

	if(rows == 0)
		rows = 1;

	if(columns == 0)
		columns = 1;

	FLOATTYPE **OldMatrix = mMatrix;

	mMatrix = new FLOATTYPE *[rows];

	assert_t(mMatrix && "Error: not enough memory!");

	if(!mMatrix)
	{
		mMatrix = OldMatrix;

		return false;
	}

	for(unsigned long i = 0; i < rows; i++)
	{
		mMatrix[i] = new FLOATTYPE[columns];

		assert_t(mMatrix[i] && "Error: not enough memory!");

		if(!mMatrix[i])
		{
			mMatrix = OldMatrix;

			return false;
		}
	}

	//Free old memory.

	for(unsigned long i = 0; i < mRows; i++)
		delete [] OldMatrix[i];

	delete [] OldMatrix;

	//Save new dims.

	mRows = rows;
	mColumns = columns;

	return true;
}

void CMatrixMN::GetDims(unsigned long &rows, unsigned long &columns) const
{
	rows = mRows;
	columns = mColumns;
}

unsigned long CMatrixMN::GetRows(void) const
{
	return mRows;
}

unsigned long CMatrixMN::GetColumns(void) const
{
	return mColumns;
}

void CMatrixMN::SwapRows(unsigned long i1, unsigned long i2)
{
	assert_t(i1 < mRows && i2 < mRows && "Index range check error.");

	if(i1 == i2)
		return;

	for(unsigned long j = 0; j < mColumns; j++)
	{
		FLOATTYPE t       = mMatrix[i1][j];
		mMatrix[i1][j] = mMatrix[i2][j];
		mMatrix[i2][j] = t;
	}
}
