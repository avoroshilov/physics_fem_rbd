#ifndef _M_VECTOR3_H_
#define _M_VECTOR3_H_

#include "all.h"

class Vec3
{
public:

	//Предоставляем на выбор два способа обращения к координатам вектора:
	union
	{
		struct
		{
			FLOATTYPE x, y, z;
		};

		struct
		{
			FLOATTYPE v[3];
		};
	};

	Vec3(void){}										//Не инициализировать. WARNING!!!
	Vec3(FLOATTYPE ix, FLOATTYPE iy, FLOATTYPE iz);				//Инициализировать по координатам.
	//~Vec3(void){}										//Деструктор. Ничего не делает.

	void				Assign(FLOATTYPE ix, FLOATTYPE iy, FLOATTYPE iz);	//Присвоить координаты (to assign a value to the variable == приписать значение переменной).
	void				Set(FLOATTYPE ix, FLOATTYPE iy, FLOATTYPE iz);

	FLOATTYPE			GetLength(void) const;						//Получение длины вектора.
	FLOATTYPE			GetSquareLength(void) const;				//Получение квадрата длины вектора.

	void				Normalize(void);							//Нормализовать вектор.
	Vec3				GetNormalized(void) const;					//Получение нормализованного вектора без изменения исходного.
	void				MakeInverse(void);							//Инвертировать вектор.

	FLOATTYPE			Dot(const Vec3 &arg) const;				//Точечное произведение векторов.
	Vec3				Cross(const Vec3 &arg) const;			//Векторное произведение векторов.
	FLOATTYPE			Triple(const Vec3 &b, const Vec3 &c) const;//Смешанное произведение векторов
	Vec3				ProjectOnto(const Vec3 &arg) const;		//Проецирование вектора на ось.

	void				TangentSpace(Vec3 &Tan1, Vec3 &Tan2) const;

	FLOATTYPE &			operator[](unsigned long index);			//Доступ к координатам по индексу.
	const FLOATTYPE &	operator[](unsigned long index) const;		//То же, но для константного экземпляра класса.

	const Vec3 &		operator =(const Vec3 &arg);			//Операция присваивания.

	bool				operator ==(const Vec3 &arg) const;		//Операция сравнения.
	bool				operator !=(const Vec3 &arg) const;		//Операция сравнения.

	bool				IsEqualAt(const Vec3 &arg, const FLOATTYPE &precision) const;

	bool				operator <(const Vec3 &arg) const;		//Defined for sorting vectors.

	Vec3				operator -(void)const;						//Унарный минус.

	Vec3				operator +(const Vec3 &arg) const;		//Оператор сложения с вектором.
	Vec3				operator -(const Vec3 &arg) const;		//Оператор вычитания вектора.
	Vec3				operator *(const FLOATTYPE &arg) const;		//Оператор умножения на число.
	Vec3				operator /(const FLOATTYPE &arg) const;		//Оператор деления на число.

	Vec3 &				operator +=(const Vec3 &arg);			//Добавление вектора.
	Vec3 &				operator -=(const Vec3 &arg);			//Вычитание вектора.
	Vec3 &				operator *=(const FLOATTYPE &arg);				//Умножение на число.
	Vec3 &				operator /=(const FLOATTYPE &arg);				//Деление на число.
};

Vec3 operator *(const FLOATTYPE &scalar, const Vec3 &vector);

const Vec3 NullVec3(0.0f, 0.0f, 0.0f);

#ifdef INLINE_BUILD
	#include "m_Vector3.inl"
#endif

#endif
