#ifndef _M_MATRIX3_H_
#define _M_MATRIX3_H_

#include "all.h"

#include "m_Vector3.h"

class CMatrix3
{
public:
	union
	{
		struct
		{
			FLOATTYPE mMatrix[3][3];			//Первый индекс определяет строку, второй - столбец.
		};

		struct
		{
			Vec3 mRows[3];
		};
	};

	CMatrix3();																// Empty default c-tor
	CMatrix3(const CMatrix3 &rMatrix);										//Конструктор копирования. Нах тут не нужен.
	CMatrix3(const Vec3 &Axis, float Angle);								// Constructs rotational matrix from axis-angle.
	//~CMatrix3(void){}														//Деструктор. Ничего не делает.

	void									MakeZero(void);								//Обнулить матрицу.
	void									MakeIdentity(void);							//Сделать матрицу матрицей идентичности.
	void									Orthogonalize(void);						//Сделать матрицу ортогональной.

	void									Transpose(void);
	CMatrix3								GetTransposed(void) const;

	CMatrix3								GetInverse(FLOATTYPE *det = NULL);

	void									ToColumnMatrix4(float m[16]) const;				//Преобразовать к формату матрицы в OpenGL.
	
	CMatrix3 								Sqrt(unsigned int niter = 10) const;
	void									PolarDecompose(CMatrix3 &p, CMatrix3 &u, unsigned int niter = 10) const;

	FLOATTYPE								Det() const;

	FLOATTYPE *								operator[](unsigned long index);			//Доступ к указателю на первый элемент строки матрицы по индексу.
	const FLOATTYPE *							operator[](unsigned long index) const;		//Доступ к указателю на первый элемент строки матрицы по индексу для константного экземпляра класса.

	CMatrix3 &								operator =(const CMatrix3 &arg);			//Оператор присваивания.

	CMatrix3								operator +(const CMatrix3 &arg) const;		//Оператор сложения матриц.
	CMatrix3								operator -(const CMatrix3 &arg) const;		//Оператор вычитания матриц.
	CMatrix3								operator *(FLOATTYPE arg) const;				//Оператор умножения матрицы на число.
	Vec3								operator *(const Vec3 &arg) const;		//Оператор умножения матрицы на вектор (поворачивает вектор).
	CMatrix3								operator *(const CMatrix3 &arg) const;		//Оператор умножения матриц.

	void									operator +=(const CMatrix3 &arg);			//Оператор добавления матрицы.
	void									operator -=(const CMatrix3 &arg);			//Оператор вычитания матрицы.
	void									operator *=(FLOATTYPE arg);					//Оператор умножения на чтсло.
	void									operator *=(const CMatrix3 &arg);			//Оператор умножения на матрицу.
};

CMatrix3 operator *(FLOATTYPE scalar, const CMatrix3 &matrix);
Vec3 operator *(const Vec3 &vec, const CMatrix3 &mat);

#ifdef INLINE_BUILD
	#include "m_Matrix3.inl"
#endif

#endif
