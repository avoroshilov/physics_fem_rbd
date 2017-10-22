#ifndef _M_VECTOR3_H_
#define _M_VECTOR3_H_

#include "all.h"

class Vec3
{
public:

	//������������� �� ����� ��� ������� ��������� � ����������� �������:
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

	Vec3(void){}										//�� ����������������. WARNING!!!
	Vec3(FLOATTYPE ix, FLOATTYPE iy, FLOATTYPE iz);				//���������������� �� �����������.
	//~Vec3(void){}										//����������. ������ �� ������.

	void				Assign(FLOATTYPE ix, FLOATTYPE iy, FLOATTYPE iz);	//��������� ���������� (to assign a value to the variable == ��������� �������� ����������).
	void				Set(FLOATTYPE ix, FLOATTYPE iy, FLOATTYPE iz);

	FLOATTYPE			GetLength(void) const;						//��������� ����� �������.
	FLOATTYPE			GetSquareLength(void) const;				//��������� �������� ����� �������.

	void				Normalize(void);							//������������� ������.
	Vec3				GetNormalized(void) const;					//��������� ���������������� ������� ��� ��������� ���������.
	void				MakeInverse(void);							//������������� ������.

	FLOATTYPE			Dot(const Vec3 &arg) const;				//�������� ������������ ��������.
	Vec3				Cross(const Vec3 &arg) const;			//��������� ������������ ��������.
	FLOATTYPE			Triple(const Vec3 &b, const Vec3 &c) const;//��������� ������������ ��������
	Vec3				ProjectOnto(const Vec3 &arg) const;		//������������� ������� �� ���.

	void				TangentSpace(Vec3 &Tan1, Vec3 &Tan2) const;

	FLOATTYPE &			operator[](unsigned long index);			//������ � ����������� �� �������.
	const FLOATTYPE &	operator[](unsigned long index) const;		//�� ��, �� ��� ������������ ���������� ������.

	const Vec3 &		operator =(const Vec3 &arg);			//�������� ������������.

	bool				operator ==(const Vec3 &arg) const;		//�������� ���������.
	bool				operator !=(const Vec3 &arg) const;		//�������� ���������.

	bool				IsEqualAt(const Vec3 &arg, const FLOATTYPE &precision) const;

	bool				operator <(const Vec3 &arg) const;		//Defined for sorting vectors.

	Vec3				operator -(void)const;						//������� �����.

	Vec3				operator +(const Vec3 &arg) const;		//�������� �������� � ��������.
	Vec3				operator -(const Vec3 &arg) const;		//�������� ��������� �������.
	Vec3				operator *(const FLOATTYPE &arg) const;		//�������� ��������� �� �����.
	Vec3				operator /(const FLOATTYPE &arg) const;		//�������� ������� �� �����.

	Vec3 &				operator +=(const Vec3 &arg);			//���������� �������.
	Vec3 &				operator -=(const Vec3 &arg);			//��������� �������.
	Vec3 &				operator *=(const FLOATTYPE &arg);				//��������� �� �����.
	Vec3 &				operator /=(const FLOATTYPE &arg);				//������� �� �����.
};

Vec3 operator *(const FLOATTYPE &scalar, const Vec3 &vector);

const Vec3 NullVec3(0.0f, 0.0f, 0.0f);

#ifdef INLINE_BUILD
	#include "m_Vector3.inl"
#endif

#endif
