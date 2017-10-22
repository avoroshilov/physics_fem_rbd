#include "m_MatrixN.h"

#pragma message(__FILE__": warning: Create copy constructors for matrices!!!")

using namespace math;

CMatrixN::CMatrixN(void)
{
	mMatrix = new FLOATTYPE *[1];
	mMatrix[0] = new FLOATTYPE[1];

	mSize = 1;
}

CMatrixN::CMatrixN(unsigned long size)
{
	assert(size > 0);

	if(size == 0)
		size = 1;

	mMatrix = new FLOATTYPE *[size];

	assert_t(mMatrix && "Error: not enough memory!");

	for(unsigned long i = 0; i < size; i++)
	{
		mMatrix[i] = new FLOATTYPE[size];
		assert_t(mMatrix[i]);
	}

	mSize = size;
}

CMatrixN::CMatrixN(const CMatrixN &arg)
{
	mSize = 0;
	mMatrix = NULL;
	*this = arg;
}

CMatrixN::CMatrixN(const CMatrix3 &arg)
{
	mSize = 0;
	mMatrix = NULL;
	*this = arg;
}

CMatrixN::CMatrixN(const CMatrix4 &arg)
{
	mSize = 0;
	mMatrix = NULL;
	*this = arg;
}


CMatrixN::~CMatrixN(void)
{
	if(mMatrix)
	{
		for(unsigned long i = 0; i < mSize; i++)
		{
			assert(mMatrix[i]);
			delete [] mMatrix[i];
		}

		delete [] mMatrix;

		mMatrix = NULL;
	}
}

/*FLOATTYPE &CMatrixN::operator ()(unsigned long i, unsigned long j)
{
//	assert_t(i < mSize && j < mSize && "Index range check error.");

	return mMatrix[i][j];
}

const FLOATTYPE &CMatrixN::operator ()(unsigned long i, unsigned long j) const
{
//	assert_t(i < mSize && j < mSize && "Index range check error.");

	return mMatrix[i][j];
}*/

CMatrixN &CMatrixN::operator -(void)
{
	for(unsigned long i = 0; i < mSize; i++)
	{
		for(unsigned long j = 0; j < mSize; j++)
		{
			mMatrix[i][j] *= -1.0;
		}
	}

	return *this;
};

CVectorN CMatrixN::operator *(const CVectorN &arg) const
{
	assert_t(arg.GetSize() == mSize);

	CVectorN result(mSize);

	if(arg.GetSize() != mSize)
	{
	//	printf("Size check failed in CMatrixN.\r\n");
		assert(arg.GetSize() == mSize && "Size mismatch!");

		return result;
	}

	for(unsigned long i = 0; i < mSize; i++)
	{
		result[i] = 0.0f;

		for(unsigned long j = 0; j < mSize; j++)
		{
			result[i] += mMatrix[i][j] * arg[j];
		}
	}

	return result;
}

CMatrixN& CMatrixN::operator =(const CMatrixN &arg)
{
//	FILE *pFile = fopen("matrix.txt","a");
//	fprintf(pFile,"\n");

//	message("a");

	if(mSize != arg.GetSize())
		this->SetSize(arg.GetSize());

	for(unsigned long i = 0; i < mSize; i++)
		for(unsigned long j = 0; j < mSize; j++)
		{
			mMatrix[i][j] = arg.mMatrix[i][j];
//		fprintf(pFile, "\t %f", mMatrix[i][j]);
		}
//	fclose(pFile);

	return *this;
}

CMatrixN & CMatrixN::operator =(const CMatrix3 &arg)
{
	if(mSize != 3)
		this->SetSize(3);

	for(unsigned long i = 0; i < mSize; i++)
		for(unsigned long j = 0; j < mSize; j++)
			mMatrix[i][j] = arg[i][j];


	return *this;
}

CMatrixN & CMatrixN::operator =(const CMatrix4 &arg)
{
	if(mSize != 4)
		this->SetSize(4);

	for(unsigned long i = 0; i < mSize; i++)
		for(unsigned long j = 0; j < mSize; j++)
			mMatrix[i][j] = arg[i][j];


	return *this;
}

CMatrixN::operator CMatrix3()
{
	CMatrix3 res;

	if(mSize != 3)
	{
		assert(false && "Dimentions mismatch in CMatrix3 CMatrixN::operator CMatrix3()");
		CMatrix3 dummy;
		dummy.MakeZero();

		return dummy;
	}

	for(unsigned long i = 0; i < 3; i++)
		for(unsigned long j = 0; j <3;j++)
		{
			res[i][j] = mMatrix[i][j];
		}

	return res;
}

CMatrixN::operator CMatrix4()
{
	CMatrix4 res;

	if(mSize != 4)
	{
		assert(false && "Dimentions mismatch in CMatrix4 CMatrixN::operator CMatrix4()");
		CMatrix4 dummy;
		dummy.MakeZero();

		return dummy;
	}

	for(unsigned long i = 0; i < 4; i++)
		for(unsigned long j = 0; j < 4;j++)
		{
			res[i][j]= mMatrix[i][j];
		}

	return res;
}


CMatrixN CMatrixN::operator+(const CMatrixN &arg) const
{
	assert_t((arg.mSize == mSize) && "Index range check error.");
	assert_t(mMatrix);
	CMatrixN result(mSize);
	for(unsigned long i = 0;i < mSize;i++)
		for(unsigned long j = 0;j < mSize;j++)
			result.mMatrix[i][j] = arg.mMatrix[i][j] + mMatrix[i][j];

	return result;
}

CMatrixN CMatrixN::operator-(const CMatrixN &arg) const
{
	assert_t((arg.mSize == mSize) && "Index range check error.");
	assert_t(mMatrix);
	CMatrixN result(mSize);
	for(unsigned long i = 0;i < mSize;i++)
		for(unsigned long j = 0;j < mSize;j++)
			result.mMatrix[i][j] = mMatrix[i][j] - arg.mMatrix[i][j];

	return result;
}

CMatrixN CMatrixN::operator *(FLOATTYPE arg) const
{
	CMatrixN result(mSize);

	for(unsigned long i = 0; i < mSize; i++)
		for(unsigned long j = 0; j < mSize; j++)
			result.mMatrix[i][j] = mMatrix[i][j] * arg;

	return result;
}

CMatrixN CMatrixN::operator *(const CMatrixN &arg) const
{
	CMatrixN result(mSize);
	result.MakeZero();
	for (unsigned long i = 0; i < mSize; i++)
		for (unsigned long j = 0; j < mSize; j++)
			for (unsigned long k = 0; k < mSize; k++)
			{
				result.mMatrix[i][j] += (mMatrix[i][k] * arg.mMatrix[k][j]);
			}

	return result;
}

CMatrixN CMatrixN::operator +=(const CMatrixN &arg)
{
	assert_t((arg.mSize == mSize) && "Index range check error.");
	assert_t(mMatrix);
	for(unsigned long i = 0;i < mSize; i++)
		for(unsigned long j = 0;j < mSize; j++)
			mMatrix[i][j] += arg.mMatrix[i][j];

	return *this;
}

CMatrixN CMatrixN::operator -=(const CMatrixN &arg)
{
	assert_t((arg.mSize == mSize) && "Index range check error.");
	assert_t(mMatrix);
	for(unsigned long i = 0; i < mSize; i++)
		for(unsigned long j = 0;j < mSize;j++)
			mMatrix[i][j] -= arg.mMatrix[i][j];

	return *this;
}

CMatrixN CMatrixN::operator *=(FLOATTYPE arg)
{

	for(unsigned long i = 0; i < mSize; i++)
		for(unsigned long j = 0; j < mSize; j++)
		mMatrix[i][j] *= arg;

	return *this;
}

CMatrixN CMatrixN::operator *=(const CMatrixN &arg)
{
	CMatrixN result(mSize);
	result.MakeZero();
	for (unsigned long i = 0; i < mSize; i++)
		for (unsigned long j = 0; j < mSize; j++)
			for (unsigned long k = 0; k < mSize; k++)
			{
				result.mMatrix[i][j] += (mMatrix[i][k] * arg.mMatrix[k][j]);
			}
	*this = result;

	return result;
}


bool CMatrixN::SetSize(unsigned long size)
{
	assert(size > 0);

	if(size == 0)
		size = 1;

	FLOATTYPE **OldMatrix = mMatrix;

	mMatrix = new FLOATTYPE *[size];

	assert_t(mMatrix && "Error: not enough memory!");

	if(!mMatrix)
	{
		mMatrix = OldMatrix;

		return false;
	}

	for(unsigned long i = 0; i < size; i++)
	{
		mMatrix[i] = new FLOATTYPE[size];

		assert_t(mMatrix[i] && "Error: not enough memory!");

		if(!mMatrix[i])
		{
			mMatrix = OldMatrix;

			return false;
		}
	}

	//Free old memory.

	for(unsigned long i = 0; i < mSize; i++)
		delete [] OldMatrix[i];

	delete [] OldMatrix;

	//Save new dims.

	mSize = size;

	return true;
}

unsigned long CMatrixN::GetSize(void) const
{
	return mSize;
}

void CMatrixN::MakeZero(void)
{
	for(unsigned long i = 0; i < mSize; i++)
		for(unsigned long j = 0; j < mSize; j++)
			mMatrix[i][j] = 0;
	
	return;
}

void CMatrixN::MakeIdentity(void)
{
	for(unsigned long i = 0; i < mSize; i++) for(unsigned long j = 0; j < mSize; j++)
	{
		if(i == j)
			mMatrix[i][j] = 1;
		else
			mMatrix[i][j] = 0;
	}

	return;
}
/*********************************************************************************************************************\
 Algorithm uses Gram-Schmidt orthogonalization.  If 'this' matrix is M = [m0|m1|m2], then orthonormal
 output matrix is Q = [q0|q1|q2],

  q0 = m0/|m0|
  q1 = (m1-(q0*m1)q0)/|m1-(q0*m1)q0|
  q2 = (m2-(q0*m2)q0-(q1*m2)q1)/|m2-(q0*m2)q0-(q1*m2)q1|

 where |V| indicates length of vector V and A*B indicates dot product of vectors A and B.

  Думаю, что базис из столбцов матрицы получится ортонормированный. TODO: проверить.
\*********************************************************************************************************************/


void CMatrixN::Orthogonalize(void)
{
	CVectorN vtemp(mSize), vtemp2(mSize), vbasissum(mSize);

	for(unsigned long col = 0; col < mSize; col++)
	{

		for(unsigned long i = 0; i < mSize; i++)
			vtemp[i] = mMatrix[i][col];
		
		vbasissum.MakeZero();
		
		for(unsigned long j = 0; j < col; j++)
		{
			CVectorN v(mSize);
			v = this->GetColumn(j);
			CVectorN vv = vtemp.ProjectOnto(v);
		//	FLOATTYPE f =  vtemp.Dot(v) / vtemp.GetSquareLength();

			vbasissum = vbasissum   + vv;
		}


		vtemp2 = (vtemp-vbasissum).GetNormalized();

		for(unsigned long i = 0; i < mSize; i++)
			mMatrix[i][col] = vtemp2[i];
	}

	return;
}


void CMatrixN::Transpose(void)
{
	FLOATTYPE t;
	for(unsigned long i = 0; i < mSize;i++) for(unsigned long j = i + 1; j < mSize;j++)
	{
		t = mMatrix[i][j];
		mMatrix[i][j]= mMatrix[j][i];
		mMatrix[j][i]= t;
	}

	return;
}

CMatrixN CMatrixN::GetTransposed(void)
{
	CMatrixN res = *this;
	res.Transpose();

	return res;
}

FLOATTYPE CMatrixN::Supplement(unsigned long row, unsigned long col)
{
//	message("%d", mSize);

	CMatrixN temp(mSize-1);
	char shiftr = 0;

	for(unsigned long i = 0; i< mSize;i++)
	{
		if(i == row)
			continue;

		shiftr = (i > row)? -1 : 0;
		char shiftc = 0;

		for(unsigned long j = 0; j < mSize;j++)
		{

			if(j == col)
				continue;
			
			shiftc = (j > col)? -1 : 0;
			
			temp(i + shiftr, j + shiftc) = mMatrix[i][j];
		}
	}

	char sign = (char) ( (row+col) % 2);

	return (sign ? -1 : 1) * temp.Det();
}

FLOATTYPE CMatrixN::Det()
{
	if(mSize == 1)
		return mMatrix[0][0];

	FLOATTYPE det = 0.0;

	for(unsigned long i = 0; i < mSize; i++)
		det += (mMatrix[0][i] * this->Supplement(0, i));

	return det;
}

CMatrixN CMatrixN::GetInverse(FLOATTYPE *det)
{
	CMatrixN inv(mSize);
	inv.MakeZero();
	FLOATTYPE Determ = this->Det();

	if(det)
		*det = Determ;

	if(!M_EQUAL(Determ, 0.0f))
	{
		FLOATTYPE InvDeterm = 1/ Determ;

		for(unsigned long i = 0; i < mSize; i++)
			for(unsigned long j = 0; j < mSize; j++)
				inv(j, i) = this->Supplement(i, j) * InvDeterm;
	}

	return inv;
}

void CMatrixN::SwapRows(unsigned long i1, unsigned long i2)
{
	assert_t(i1 < mSize && i2 < mSize && "Index range check error.");

	if(i1 == i2)
		return;

	for(unsigned long j = 0; j < mSize; j++)
	{
		FLOATTYPE t        = mMatrix[i1][j];
		mMatrix[i1][j] = mMatrix[i2][j];
		mMatrix[i2][j] = t;
	}
}

bool CMatrixN::SetColumn(CVectorN column, unsigned long index)
{
	if(!(column.GetSize() == mSize))
		return false;

	assert_t(mMatrix);

	if(index >= mSize)
		return false;

	for(unsigned long i = 0; i < mSize; i++)
		mMatrix[i][index] = column[i];

	return true;
}

bool CMatrixN::SetRow(CVectorN row, unsigned long index)
{
	if(!(row.GetSize() == mSize))
		return false;

	assert_t(mMatrix);

	if(index >= mSize)
		return false;

	for(unsigned long i = 0; i < mSize; i++)
		mMatrix[index][i] = row[i];

	return true;
}

bool CMatrixN::Set(unsigned long i, unsigned long j, FLOATTYPE fvalue)
{
	assert_t(mMatrix);

	if(i >= mSize)
		return false;

	if(j >= mSize)
		return false;

	mMatrix[i][j] = fvalue;

	return true;
}


CVectorN CMatrixN::GetColumn(unsigned long index)
{
	CVectorN result(mSize);

	for(unsigned long i = 0; i<mSize; i++)
		result[i] = mMatrix[i][index];

	return result;
}

CVectorN CMatrixN::GetRow(unsigned long index)
{
	CVectorN result(mSize);

	for(unsigned long i = 0; i<mSize; i++)
		result[i] = mMatrix[index][i];

	return result;
}

void CMatrixN::Log(void)
{
	char str[1024*10];

	for(unsigned long r = 0; r < mSize; r++)
	{
		str[0] = 0;

		for(unsigned long c = 0; c < mSize; c++)
		{
			char s[200];

			sprintf(s, "%.4f ", mMatrix[r][c]);
			strcat(str, s);
		}

		gLog.Print("%s", str);
	}
}

/*CMatrixN operator*(FLOATTYPE scalar, const CMatrixN &matrix)
{
	return matrix * scalar;
}*/
