#include "m_GaussSolver.h"

#include "m_Const.h"

#include "helpers/debug.h"

#pragma message(__FILE__" : warning: Turn on assert macros!")

//

using namespace math;

bool CGaussSolver::SetSize(unsigned long size)
{
	return mMatrix.SetSize(size) && mVector.SetSize(size) && mResult.SetSize(size);
}

bool CGaussSolver::Solve(void)
{
	assert( (mMatrix.GetSize() == mVector.GetSize()) && (mMatrix.GetSize() == mResult.GetSize()) );

	if( (mMatrix.GetSize() != mVector.GetSize()) || (mMatrix.GetSize() != mResult.GetSize()) )
		return false;

//	mResult.MakeZero();
	const unsigned long size = mMatrix.GetSize();
	//������ ��������. ��������� � ���� ��������� ����� �������� ������� ������,
	//indices[����� �������] = ����� ����������, ������� ��� ���� ��������
	unsigned long* indices = new unsigned long[size];
	
	for(unsigned long i = 0; i<size;i++)
		indices[i] = i;

//	com::CArray<bool> OmittedVar;
//	OmittedVar.SetUsedSize(size);
//	com::CArray<unsigned long> StartVar;
//	StartVar.SetUsedSize(size);
//	unsigned long shift = 0;
/*	for(unsigned long i = 0; i< size; i++)
	{
		StartVar[i]= i;
		OmittedVar[i]= false;
		mResult[i]= 0.0f;
	}*/

	//-----------------------------------------------------------------------------
	// ������ ������
	//-----------------------------------------------------------------------------
	//������� ���� - ��� ������ ����������. � ���� ����������� �����������
	//1 ����������
	unsigned long var = 0;

	for(var = 0;var<size;var++)
	{
		//� ������� ����� ��� ����� ������� (�� ����������) � ���� �� �����
		//��������� ���������
		unsigned long col;
		for(col = var;col<size;col++)
		{
			//�� ���������� ����� ��� ����� ������ � ��������� ��������� � ������ ������

			unsigned long row;
			for(row = var;row<size;row++)
			{
				//���� ����� ������� ������...
				if(mMatrix(row, col)!= 0)
				{
					//� ���� �� �� ������ ��������� (�.�. ������� � ������,
					//��������������� ������ ����������� ����������), �� ����
					//�������� ������� �������
					if(row!= var)
					{
						FLOATTYPE temp;

						for(unsigned long i = var;i<size;i++)
						{
							temp = mMatrix(var, i);
							mMatrix(var, i)= mMatrix(row, i);
							mMatrix(row, i) = temp;
						}

						temp = mVector[var];
						mVector[var]= mVector[row];
						mVector[row]= temp;
					}
				//	���� ������ ��������� ������� (� ������ � ��� ������������
				// � ������ �����, �� ������ ������������� �������� ������ �����
					//������� ������� �������

					break;
				}
			}
			//������������ ����� ������� ��������, ��� � ������ ������� ���
			//������ ��������� �������
			if(row!= size)
			{
				//���� ���� ������� - �� ������ ���������, �.�. ���� ��������� �������
				//�������, ����� ��� �������� ��������� �� ������ �����, �����
				//���������������� ������� ���� �� ����� �����������
				if(col!= var)
				{
					FLOATTYPE temp;
					for(unsigned long i = 0; i<size;i++)
					{
						temp = mMatrix(i, var);
						mMatrix(i, var) = mMatrix(i, col);
						mMatrix(i, col) = temp;
					}
					//�� ���� ����� ������� ������������� ����� ����� ����������� ��������
					//� �����������, �� ������� �� ��������.
					unsigned long ltemp = indices[var];
					indices[var] = indices[col];
					indices[col]= ltemp;
				}
				//���� ������ ������� � ��������� ���������, ������ �� ����
				break;
			}
		}
		//��� ������� ��������, ��� ������� � ��������� ��������� �� ������,
		//�� ���� �� ����� ����������� �������, �� ��� ������� ������� � ���� �������
		//��������� �������, �.�. �� ����� ���� ���������, ���� ���� �������.
		//� ����� ������ ��� ������� ���������� ������� ������. ��������� ������� ����
		if(col == size)
		{
			break;
		}
		//���� �� �����, �� �� ������� (var, var) ����� ��������� �������
		//��������� ��� ���� ��������� �������.
		//FLOATTYPE div = 1/mMatrix(var, var);
		//��� 1. ������� ���� ����� ��� ����������� �������
		FLOATTYPE d = mMatrix(var, var);
		mVector[var] = mVector[var] / d;

		for(unsigned long i = var; i<size;i++)
			mMatrix(var, i) = mMatrix(var, i) / d;

		//��� 2. ���������� ������ ������ ( � �������������) �� ���� ����� ���� ������.
		for(unsigned long i = var+1; i<size;i++)
			mVector[i] -= mVector[var] * mMatrix(i, var);

		for(unsigned long j = var+1; j<size;j++)
		{
			FLOATTYPE sub = mMatrix(j, var);
			for(unsigned long i = var; i<size;i++)
				mMatrix(j, i) -= mMatrix(var, i)*sub;
		}
	}
	//*******************************
	//�������� ������� �������
	//*******************************

	//��� ������� ������, ��� ����� ����������� ��������� ������ ������������
	//����� ����������. �.�. ����� ������� ���� ������� ������. ���� ��������
	//�������, ��������������� ���� ������� ���������, �� ������� ���, �����
	//������� ���������.
	if(var != size)
	{
		for(unsigned long i = var; i<size; i++)
		{
			//���������� ����� ����������, ������ ��� "������ � ��� �������"
			//�� ��� �������, ��� ���������� ������� ������
			if(!M_EQUAL_AT(mVector[i], 0, 0.000001f))
			{
				delete[] indices;

				return false;
			}
		}
	}

	//******************************
	//�������� ��� ������
	//******************************

	//����� ��� ������ ��� - ������ ����������, ��� ������� ������ ����� �������
	//��������� �� ��������� ���������� �� ����.
	for(unsigned long i = var;i>0;i--)
	{
	//	v[i] /= m[i*dim + i];
		//��������� ������� ��� ������ ���������
		//-1 ����� ������������, ������ ��� ���������� ����� ������ 1 ������
		//"����������� ������", ��� ����, ���� ���� ������� �� ������ � ����� -
		//���� ������������� ������.
		for(unsigned long j = i-1; j>0;j--)
		{

			mVector[j-1] -= mMatrix(j-1, i-1)*mVector[i-1];
		}
	}
	//**************************************
	//������������� ������� ����������
	//**************************************
	for(unsigned long i = 0; i<size;i++)
	{
		mResult[indices[i]]= mVector[i];
	}

	delete[] indices;

	return true;
}

bool CGaussSolver::Solve0(void)
{
	assert( (mMatrix.GetSize() == mVector.GetSize()) && (mMatrix.GetSize() == mResult.GetSize()) );

	if( (mMatrix.GetSize() != mVector.GetSize()) || (mMatrix.GetSize() != mResult.GetSize()) )
		return false;
//	mResult.MakeZero();
	const unsigned long size = mMatrix.GetSize();
	//������ ��������. ��������� � ���� ��������� ����� �������� ������� ������,
	//indices[����� �������] = ����� ����������, ������� ��� ���� ��������
	unsigned long* indices = new unsigned long[size];
	
	for(unsigned long i = 0; i<size;i++)
		indices[i] = i;

//	com::CArray<bool> OmittedVar;
//	OmittedVar.SetUsedSize(size);
//	com::CArray<unsigned long> StartVar;
//	StartVar.SetUsedSize(size);
//	unsigned long shift = 0;
/*	for(unsigned long i = 0; i< size; i++)
	{
		StartVar[i]= i;
		OmittedVar[i]= false;
		mResult[i]= 0.0f;
	}*/

	//-----------------------------------------------------------------------------
	// ������ ������
	//-----------------------------------------------------------------------------
	//������� ���� - ��� ������ ����������. � ���� ����������� �����������
	//1 ����������
	unsigned long var = 0;
	for(var = 0;var<size;var++)
	{
		//� ������� ����� ��� ����� ������� (�� ����������) � ���� �� �����
		//��������� ���������
		unsigned long col;
		for(col = var;col<size;col++)
		{
			//�� ���������� ����� ��� ����� ������ � ��������� ��������� � ������ ������

			unsigned long row;
			for(row = var;row<size;row++)
			{
				//���� ����� ������� ������...
				if(mMatrix(row, col)!= 0)
				{
					//� ���� �� �� ������ ��������� (�.�. ������� � ������,
					//��������������� ������ ����������� ����������), �� ����
					//�������� ������� �������
					if(row!= var)
					{
						FLOATTYPE temp;
						for(unsigned long i = var;i<size;i++)
						{
							temp = mMatrix(var, i);
							mMatrix(var, i)= mMatrix(row, i);
							mMatrix(row, i) = temp;
						}
						temp = mVector[var];
						mVector[var]= mVector[row];
						mVector[row]= temp;
					}
				//	���� ������ ��������� ������� (� ������ � ��� ������������
				// � ������ �����, �� ������ ������������� �������� ������ �����
					//������� ������� �������

					break;
				}
			}
			//������������ ����� ������� ��������, ��� � ������ ������� ���
			//������ ��������� �������
			if(row!= size)
			{
				//���� ���� ������� - �� ������ ���������, �.�. ���� ��������� �������
				//�������, ����� ��� �������� ��������� �� ������ �����, �����
				//���������������� ������� ���� �� ����� �����������
				if(col!= var)
				{
					FLOATTYPE temp;
					for(unsigned long i = 0; i<size;i++)
					{
						temp = mMatrix(i, var);
						mMatrix(i, var) = mMatrix(i, col);
						mMatrix(i, col) = temp;
					}
					//�� ���� ����� ������� ������������� ����� ����� ����������� ��������
					//� �����������, �� ������� �� ��������.
					unsigned long ltemp = indices[var];
					indices[var] = indices[col];
					indices[col]= ltemp;
				}
				//���� ������ ������� � ��������� ���������, ������ �� ����
				break;
			}
		}
		//��� ������� ��������, ��� ������� � ��������� ��������� �� ������,
		//�� ���� �� ����� ����������� �������, �� ��� ������� ������� � ���� �������
		//��������� �������, �.�. �� ����� ���� ���������, ���� ���� �������.
		//� ����� ������ ��� ������� ���������� ������� ������. ��������� ������� ����
		if(col == size)
		{
			break;
		}
		//���� �� �����, �� �� ������� (var, var) ����� ��������� �������
		//��������� ��� ���� ��������� �������.
		FLOATTYPE div = 1/mMatrix(var, var);
		//��� 1. ������� ���� ����� ��� ����������� �������
		mVector[var] *= div;

		for(unsigned long i = var; i<size;i++)
			mMatrix(var, i) *= div;

		//��� 2. ���������� ������ ������ ( � �������������) �� ���� ����� ���� ������.
		for(unsigned long i = var+1; i<size;i++)
			mVector[i] -= mVector[var] * mMatrix(i, var);

		for(unsigned long j = var+1; j<size;j++)
		{
			FLOATTYPE sub = mMatrix(j, var);
			for(unsigned long i = var; i<size;i++)
				mMatrix(j, i) -= mMatrix(var, i)*sub;
		}
	}
	//*******************************
	//�������� ������� �������
	//*******************************

	//��� ������� ������, ��� ����� ����������� ��������� ������ ������������
	//����� ����������. �.�. ����� ������� ���� ������� ������. ���� ��������
	//�������, ��������������� ���� ������� ���������, �� ������� ���, �����
	//������� ���������.
	if(var != size)
	{
		for(unsigned long i = var; i<size; i++)
		{
			//���������� ����� ����������, ������ ��� "������ � ��� �������"
			//�� ��� �������, ��� ���������� ������� ������
			if(!M_EQUAL_AT(mVector[i], 0, 0.000001f))
			{
				delete[] indices;

				return false;
			}
		}
	}

	//******************************
	//�������� ��� ������
	//******************************

	//����� ��� ������ ��� - ������ ����������, ��� ������� ������ ����� �������
	//��������� �� ��������� ���������� �� ����.
	for(unsigned long i = var;i>0;i--)
	{
	//	v[i] /= m[i*dim + i];
		//��������� ������� ��� ������ ���������
		//-1 ����� ������������, ������ ��� ���������� ����� ������ 1 ������
		//"����������� ������", ��� ����, ���� ���� ������� �� ������ � ����� -
		//���� ������������� ������.
		for(unsigned long j = i-1; j>0;j--)
		{

			mVector[j-1] -= mMatrix(j-1, i-1)*mVector[i-1];
		}
	}
	//**************************************
	//������������� ������� ����������
	//**************************************
	for(unsigned long i = 0; i<size;i++)
	{
		mResult[indices[i]]= mVector[i];
	}

	delete[] indices;

	return true;
}

bool CGaussSolver::Solve1(void)
{
	assert( (mMatrix.GetSize() == mVector.GetSize()) && (mMatrix.GetSize() == mResult.GetSize()) );

	if( (mMatrix.GetSize() != mVector.GetSize()) || (mMatrix.GetSize() != mResult.GetSize()) )
		return false;

	const unsigned long size = mMatrix.GetSize();

	//-----------------------------------------------------------------------------
	// Prepare the matrix
	//-----------------------------------------------------------------------------

	for(unsigned long n = 0; n < size - 1; n++)
	{
		//-----------------------------------------------------------------------------
		// Find a base row
		//-----------------------------------------------------------------------------

		for(unsigned long i = n; i < size; i++)
		{
			//TODO: replace here with abs.
			if(mMatrix(i, n) < -0.000001f || mMatrix(i, n) > 0.000001f)
			{
				mMatrix.SwapRows(n, i);
				mVector.SwapElements(n, i);

				break;
			}

			if(i == size - 1)//Bad matrix.
				return false;
		}

		//-----------------------------------------------------------------------------
		// A pass
		//-----------------------------------------------------------------------------

		for(unsigned long i = n + 1; i < size; i++)
		{
			if(mMatrix(n, n) > -0.000001f && mMatrix(n, n) < 0.000001f)
				return false;

			FLOATTYPE k = mMatrix(i, n) / mMatrix(n, n);

			for(unsigned long j = n + 1; j < size; j++)
				mMatrix(i, j) -= mMatrix(n, j) * k;

			mVector[i] -= mVector[n] * k;
		}
	}

	//-----------------------------------------------------------------------------
	// Extract results
	//-----------------------------------------------------------------------------

	for(unsigned long i = size - 1; ; i--)
	{
		FLOATTYPE res;

		if(mMatrix(i, i) > -0.0001f && mMatrix(i, i) < 0.0001f)
		{
			mResult[i] = 1.0f;
		}
		else
		{
			res = mVector[i];

			for(unsigned long j = i + 1; j < size; j++)
				res -= mMatrix(i, j) * mResult[j];

			mResult[i] = res / mMatrix(i, i);
		}

		if(i == 0)
			break;
	}

	return true;
}