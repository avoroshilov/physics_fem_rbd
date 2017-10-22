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

	//������������� �� ����� ��� ������� ��������� � ����������� �����������:
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
	CQuaternion(void){}											//�� ����������������.
	CQuaternion(FLOATTYPE iw, FLOATTYPE ix, FLOATTYPE iy, FLOATTYPE iz);	//���������������� �� �����������.
	CQuaternion(const CQuaternion &quat);						//���������������� �� ����������� (copy constructor).
	CQuaternion(const CMatrix3 &matrix);						//���������������� �� �������.
	CQuaternion(const Vec3 &axis, FLOATTYPE angle);			//���������������� �� ���� ����-��� (������ ��� ������ ���� ������������).
	//~CQuaternion(void){}										//����������. ������ �� ������.

	void					Assign(FLOATTYPE iw, FLOATTYPE ix, FLOATTYPE iy, FLOATTYPE iz);	//Assign values to quaternion's coordinates.

	void					MakeZero(void);									//�������� ����������.
	void					MakeIdentity(void);								//������� ���������� ���������� (������� �������).
	void					Normalize(void);								//������������� ���������� (������� ���������).
	void					MakeInverse(void);								//������������� (��������� � ���������� �����������!).
	void					Conjugate(void);								//������� ���������� ����������.

	FLOATTYPE					Dot(const CQuaternion &arg) const;				//�������� ������������ ������������ (�� ���� �����).

	//���������
	FLOATTYPE &				operator[](int index);							//������ � ����������� �� �������.
	const FLOATTYPE &			operator[](int index) const;					//������ � ����������� �� ������� ��� ������������ ���������� ������.

	CQuaternion &			operator =(const CQuaternion &arg);				//�������� ������������.

	bool					operator ==(const CQuaternion &arg) const;		//�������� ���������.
	bool					operator !=(const CQuaternion &arg) const;		//�������� ���������.

	CQuaternion				operator -(void) const;							//������� �����.

	CQuaternion				operator +(const CQuaternion &arg) const;		//�������� �������� ������������.
	CQuaternion				operator -(const CQuaternion &arg) const;		//�������� ��������� ������������.
	CQuaternion				operator *(const CQuaternion &arg) const;		//�������� ��������� ������������.
	CQuaternion				operator *(FLOATTYPE fScalar) const;				//�������� ��������� ����������� �� �����.
	CQuaternion				operator /(FLOATTYPE fScalar) const;				//�������� ������� ����������� �� �����.

	Vec3				operator *(const Vec3 &vector) const;		//�������� ��������� �� ������ (������������ ������).

	CQuaternion &			operator +=(const CQuaternion &arg);			//�������� ���������� �����������.
	CQuaternion &			operator -=(const CQuaternion &arg);			//�������� ��������� �����������.
	CQuaternion &			operator *=(FLOATTYPE fScalar);					//�������� ��������� �� �����.
	CQuaternion &			operator /=(FLOATTYPE fScalar);					//�������� ������� �� �����.

	//�������������� ����� �������������, ��������� 3x3 � ������ ����-���
	void					FromMatrix3(const CMatrix3 &matrix);			//������������� ������� 3x3 � ����������.
	void					ToMatrix3(CMatrix3 &matrix) const;				//������������� ���������� � ������� 3x3.

	void					FromAxisAngle(const Vec3 &axis, FLOATTYPE angle);	//������������� ���� ����-��� � ����������.
	void					ToAxisAngle(Vec3 &axis, FLOATTYPE &angle) const;	//������������� ���������� � ���� ����-���.

	//����������� �������� ������������ - �� �����������! TODO: Test.
	static CQuaternion		Slerp(FLOATTYPE t, const CQuaternion &p, const CQuaternion &q);

	//����������� ���������� ������������ - �� �����������! TODO: Test.
	static CQuaternion		Squad(FLOATTYPE t, const CQuaternion &q0, const CQuaternion &a0, const CQuaternion &a1, const CQuaternion &q1);

protected:
	static const unsigned char mNext[3];	//���������� ��� ���������� ����������. ��������, ������ ��� ��������� ��� ���� ����������� ������.
};

CQuaternion operator *(FLOATTYPE scalar, const CQuaternion &quat);

#ifdef INLINE_BUILD
	#include "m_Quaternion.inl"
#endif

#endif
