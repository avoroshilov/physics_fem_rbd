#ifndef _M_MATRIX_N_H_
#define _M_MATRIX_N_H_
#pragma once
#include "all.h"

#include "m_VectorN.h"
#include "m_Matrix3.h"
#include "m_Matrix4.h"

namespace math
{
	class CMatrixN
	{
	public:
		CMatrixN(void);//1x1 matrix.
		CMatrixN(const CMatrixN &arg);
		CMatrixN(const CMatrix3 &arg);
		CMatrixN(const CMatrix4 &arg);

		explicit CMatrixN(unsigned long size);

		~CMatrixN(void);

		inline FLOATTYPE &						operator ()(unsigned long i, unsigned long j)
		{
			return mMatrix[i][j];
		}
		
		inline const FLOATTYPE &				operator ()(unsigned long i, unsigned long j) const
		{
			return mMatrix[i][j];
		}
		
		inline const FLOATTYPE &				Get(unsigned long i, unsigned long j) const
		{
			return mMatrix[i][j];
		}
		
		CVectorN							operator *(const CVectorN &arg) const;
		CMatrixN &							operator =(const CMatrixN &arg);			//Оператор присваивания.
		CMatrixN &							operator =(const CMatrix3 &arg);			//Оператор присваивания.
		CMatrixN &							operator =(const CMatrix4 &arg);			//Оператор присваивания.
		CMatrixN &							operator -(void);

											operator CMatrix3() ;
											operator CMatrix4() ;

		CMatrixN							operator +(const CMatrixN &arg) const;		//Оператор сложения матриц.
		CMatrixN							operator -(const CMatrixN &arg) const;		//Оператор вычитания матриц.
		CMatrixN							operator *(FLOATTYPE arg) const;				//Оператор умножения матрицы на число.
	//	CVectorN							operator *(const CVectorN &arg) const;		//Оператор умножения матрицы на вектор (поворачивает вектор).
		CMatrixN							operator *(const CMatrixN &arg) const;		//Оператор умножения матриц.

		CMatrixN							operator +=(const CMatrixN &arg);			//Оператор добавления матрицы.
		CMatrixN							operator -=(const CMatrixN &arg);			//Оператор вычитания матрицы.
		CMatrixN							operator *=(FLOATTYPE arg);						//Оператор умножения на чтсло.
		CMatrixN							operator *=(const CMatrixN &arg);			//Оператор умножения на матрицу.


		//This will not copy previous contents.
		bool								SetSize(unsigned long size);
		unsigned long						GetSize(void) const;

		bool								SetColumn(CVectorN column, unsigned long index);
		bool								SetRow(CVectorN row, unsigned long index);
		bool								Set(unsigned long i, unsigned long j, FLOATTYPE fvalue);
		CVectorN							GetColumn(unsigned long index);
		CVectorN							GetRow(unsigned long index);


		void								MakeZero(void);								//Обнулить матрицу.
		void								MakeIdentity(void);							//Сделать матрицу матрицей идентичности.
		void								Orthogonalize(void);						//Сделать матрицу ортогональной.

		void								Transpose(void);
		CMatrixN							GetTransposed(void);

		FLOATTYPE								Supplement(unsigned long row, unsigned long col);
		FLOATTYPE								Det();
		CMatrixN							GetInverse(FLOATTYPE *det = NULL);

		void								SwapRows(unsigned long i1, unsigned long i2);
		void								Log(void);
	private:
		unsigned long						mSize;
		FLOATTYPE **							mMatrix;
	};

	inline CMatrixN operator *(FLOATTYPE scalar, const CMatrixN &matrix)
	{
		return matrix * scalar;
	}
}

#endif
