#ifndef _VECTOR_N_H_
#define _VECTOR_N_H_
#pragma once
#include "all.h"
#include "helpers/debug.h"
#include "helpers/log.h"
#include "m_Const.h"
#include "m_Vector3.h"

namespace math
{
	/*
	struct MathError
	{
		char *err;
		MathError(const char * errmsg = "unspecified math error")
		{
			err = errmsg;
		}
	};
	*/
	class CVectorN
	{

	public:
		CVectorN(void);//Initialize a vector of 1 element.
		explicit CVectorN(unsigned long size);//'size' should not be 0.
		CVectorN(const CVectorN &rh);
		CVectorN(const Vec3 &rh);

		bool Set(unsigned long begin, unsigned long n,...);

		~CVectorN(void);
		FLOATTYPE								GetLength(void) const;						//Получение длины вектора.
		FLOATTYPE								GetSquareLength(void) const;				//Получение квадрата длины вектора.

		void								Normalize(void);							//Нормализовать вектор.
		CVectorN							GetNormalized(void) const;					//Получение нормализованного вектора без изменения исходного.
		void								MakeInverse(void);							//Инвертировать вектор.

		FLOATTYPE								Dot(const CVectorN &arg) const ;				//Точечное произведение векторов.
		CVectorN							ProjectOnto(const CVectorN &arg) const ;		//Проецирование вектора на ось.


		CVectorN &							operator =(const CVectorN &arg);
		CVectorN &							operator =(const Vec3 &arg);
											operator Vec3() ;
		bool								operator ==(const CVectorN &arg) const;		//Операция сравнения.
		bool								operator !=(const CVectorN &arg) const;		//Операция сравнения.

		bool								IsEqualAt(const CVectorN &arg, const FLOATTYPE &precision) const;


		inline	FLOATTYPE &					operator [](unsigned long i)
		{
			return mVector[i];
		};
		inline const FLOATTYPE &				operator [](unsigned long i) const
		{
			return mVector[i];
		};
		inline const FLOATTYPE &				Get(unsigned long i) const
		{
			return mVector[i];
		};
		CVectorN							operator -(void)const;						//Унарный минус.
		CVectorN							operator +(const CVectorN &arg) const;		//Оператор сложения с вектором.
		CVectorN							operator -(const CVectorN &arg) const;		//Оператор вычитания вектора.

		CVectorN							operator *(const FLOATTYPE &arg) const;
		CVectorN							operator /(const FLOATTYPE &arg) const;			//Оператор деления на число.

		CVectorN &							operator +=(const CVectorN &arg);
		CVectorN &							operator -=(const CVectorN &arg);			//Вычитание вектора.
		CVectorN &							operator *=(const FLOATTYPE &arg);				//Умножение на число.
		CVectorN &							operator /=(const FLOATTYPE &arg);				//Деление на число.

		void								MakeZero(void);

		//This will copy previous contents.
		bool								SetSize(unsigned long size);
		unsigned long						GetSize(void) const;

		void								SwapElements(unsigned long i, unsigned long j);

		void								Log(void);
	private:
		unsigned long						mSize;
		FLOATTYPE *							mVector;
	};
	inline CVectorN operator *(const FLOATTYPE &scalar, const CVectorN &vector)
	{
		return vector * scalar;
	}
}



#endif
