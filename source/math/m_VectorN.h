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
		FLOATTYPE								GetLength(void) const;						//��������� ����� �������.
		FLOATTYPE								GetSquareLength(void) const;				//��������� �������� ����� �������.

		void								Normalize(void);							//������������� ������.
		CVectorN							GetNormalized(void) const;					//��������� ���������������� ������� ��� ��������� ���������.
		void								MakeInverse(void);							//������������� ������.

		FLOATTYPE								Dot(const CVectorN &arg) const ;				//�������� ������������ ��������.
		CVectorN							ProjectOnto(const CVectorN &arg) const ;		//������������� ������� �� ���.


		CVectorN &							operator =(const CVectorN &arg);
		CVectorN &							operator =(const Vec3 &arg);
											operator Vec3() ;
		bool								operator ==(const CVectorN &arg) const;		//�������� ���������.
		bool								operator !=(const CVectorN &arg) const;		//�������� ���������.

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
		CVectorN							operator -(void)const;						//������� �����.
		CVectorN							operator +(const CVectorN &arg) const;		//�������� �������� � ��������.
		CVectorN							operator -(const CVectorN &arg) const;		//�������� ��������� �������.

		CVectorN							operator *(const FLOATTYPE &arg) const;
		CVectorN							operator /(const FLOATTYPE &arg) const;			//�������� ������� �� �����.

		CVectorN &							operator +=(const CVectorN &arg);
		CVectorN &							operator -=(const CVectorN &arg);			//��������� �������.
		CVectorN &							operator *=(const FLOATTYPE &arg);				//��������� �� �����.
		CVectorN &							operator /=(const FLOATTYPE &arg);				//������� �� �����.

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
