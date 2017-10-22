#ifndef _M_QUATERNION_H_
#define _M_QUATERNION_H_

#include "all.h"

#include "m_Matrix3.h"
#include "m_Vector3.h"

/*********************************************************************************************************************\
 A quaternion is q = w + x*i + y*j + z*k where(w, x, y, z) is not necessarily a unit length vector in 4D.
\*********************************************************************************************************************/
class CQuaternion
{
public:

	//Предоставляем на выбор два способа обращения к координатам кватерниона:
	union
	{
		struct
		{
			FLOATTYPE w, x, y, z;
		};

		struct
		{
			FLOATTYPE mCoords[4];
		};
	};

	//Construction and destruction
	CQuaternion(void){}											//Не инициализировать.
	CQuaternion(FLOATTYPE iw, FLOATTYPE ix, FLOATTYPE iy, FLOATTYPE iz);	//Инициализировать по координатам.
	CQuaternion(const CQuaternion &quat);						//Инициализировать по кватерниону (copy constructor).
	CQuaternion(const CMatrix3 &matrix);						//Инициализировать по матрице.
	CQuaternion(const Vec3 &axis, FLOATTYPE angle);			//Инициализировать по паре угол-ось (вектор оси должен быть нормализован).
	//~CQuaternion(void){}										//Деструктор. Ничего не делает.

	void					Assign(FLOATTYPE iw, FLOATTYPE ix, FLOATTYPE iy, FLOATTYPE iz);	//Assign values to quaternion's coordinates.

	void					MakeZero(void);									//Обнулить кватернион.
	void					MakeIdentity(void);								//Сделать кватернион идентичным (подобно матрице).
	void					Normalize(void);								//Нормализовать кватернион (сделать единичным).
	void					MakeInverse(void);								//Инвертировать (применять к единичному кватерниону!).
	void					Conjugate(void);								//Сделать кватернион сопряжённым.

	FLOATTYPE					Dot(const CQuaternion &arg) const;				//Точечное произведение кватернионов (не знаю зачем).

	//Операторы
	FLOATTYPE &				operator[](int index);							//Доступ к координатам по индексу.
	const FLOATTYPE &			operator[](int index) const;					//Доступ к координатам по индексу для константного экземпляра класса.

	CQuaternion &			operator =(const CQuaternion &arg);				//Оператор присваивания.

	bool					operator ==(const CQuaternion &arg) const;		//Оператор сравнения.
	bool					operator !=(const CQuaternion &arg) const;		//Оператор сравнения.

	CQuaternion				operator -(void) const;							//Унарный минус.

	CQuaternion				operator +(const CQuaternion &arg) const;		//Оператор сложения кватернионов.
	CQuaternion				operator -(const CQuaternion &arg) const;		//Оператор вычитания кватернионов.
	CQuaternion				operator *(const CQuaternion &arg) const;		//Оператор умножения кватернионов.
	CQuaternion				operator *(FLOATTYPE fScalar) const;				//Оператор умножения кватерниона на число.
	CQuaternion				operator /(FLOATTYPE fScalar) const;				//Оператор деления кватерниона на число.

	Vec3				operator *(const Vec3 &vector) const;		//Оператор умножения на вектор (поворачивает вектор).

	CQuaternion &			operator +=(const CQuaternion &arg);			//Оператор добавления кватерниона.
	CQuaternion &			operator -=(const CQuaternion &arg);			//Оператор вычитания кватерниона.
	CQuaternion &			operator *=(FLOATTYPE fScalar);					//Оператор умножения на число.
	CQuaternion &			operator /=(FLOATTYPE fScalar);					//Оператор деления на число.

	//Преобразования между кватернионами, матрицами 3x3 и парами угол-ось
	void					FromMatrix3(const CMatrix3 &matrix);			//Преобразовать матрицу 3x3 в кватернион.
	void					ToMatrix3(CMatrix3 &matrix) const;				//Преобразовать кватернион в матрицу 3x3.

	void					FromAxisAngle(const Vec3 &axis, FLOATTYPE angle);	//Преобразовать пару угол-ось в кватернион.
	void					ToAxisAngle(Vec3 &axis, FLOATTYPE &angle) const;	//Преобразовать кватернион в пару угол-ось.

	//Сферическая линейная интерполяция - не тестировано! TODO: Test.
	static CQuaternion		Slerp(FLOATTYPE t, const CQuaternion &p, const CQuaternion &q);

	//Сферическая кубическая интерполяция - не тестировано! TODO: Test.
	static CQuaternion		Squad(FLOATTYPE t, const CQuaternion &q0, const CQuaternion &a0, const CQuaternion &a1, const CQuaternion &q1);

protected:
	static const unsigned char mNext[3];	//Необходима для внутренних алгоритмов. Статична, потому что одинакова для всех экземпляров класса.
};

CQuaternion operator *(FLOATTYPE scalar, const CQuaternion &quat);

#ifdef INLINE_BUILD
	#include "m_Quaternion.inl"
#endif

#endif
