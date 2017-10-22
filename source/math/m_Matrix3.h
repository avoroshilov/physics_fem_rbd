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
			FLOATTYPE mMatrix[3][3];			//������ ������ ���������� ������, ������ - �������.
		};

		struct
		{
			Vec3 mRows[3];
		};
	};

	CMatrix3();																// Empty default c-tor
	CMatrix3(const CMatrix3 &rMatrix);										//����������� �����������. ��� ��� �� �����.
	CMatrix3(const Vec3 &Axis, float Angle);								// Constructs rotational matrix from axis-angle.
	//~CMatrix3(void){}														//����������. ������ �� ������.

	void									MakeZero(void);								//�������� �������.
	void									MakeIdentity(void);							//������� ������� �������� ������������.
	void									Orthogonalize(void);						//������� ������� �������������.

	void									Transpose(void);
	CMatrix3								GetTransposed(void) const;

	CMatrix3								GetInverse(FLOATTYPE *det = NULL);

	void									ToColumnMatrix4(float m[16]) const;				//������������� � ������� ������� � OpenGL.
	
	CMatrix3 								Sqrt(unsigned int niter = 10) const;
	void									PolarDecompose(CMatrix3 &p, CMatrix3 &u, unsigned int niter = 10) const;

	FLOATTYPE								Det() const;

	FLOATTYPE *								operator[](unsigned long index);			//������ � ��������� �� ������ ������� ������ ������� �� �������.
	const FLOATTYPE *							operator[](unsigned long index) const;		//������ � ��������� �� ������ ������� ������ ������� �� ������� ��� ������������ ���������� ������.

	CMatrix3 &								operator =(const CMatrix3 &arg);			//�������� ������������.

	CMatrix3								operator +(const CMatrix3 &arg) const;		//�������� �������� ������.
	CMatrix3								operator -(const CMatrix3 &arg) const;		//�������� ��������� ������.
	CMatrix3								operator *(FLOATTYPE arg) const;				//�������� ��������� ������� �� �����.
	Vec3								operator *(const Vec3 &arg) const;		//�������� ��������� ������� �� ������ (������������ ������).
	CMatrix3								operator *(const CMatrix3 &arg) const;		//�������� ��������� ������.

	void									operator +=(const CMatrix3 &arg);			//�������� ���������� �������.
	void									operator -=(const CMatrix3 &arg);			//�������� ��������� �������.
	void									operator *=(FLOATTYPE arg);					//�������� ��������� �� �����.
	void									operator *=(const CMatrix3 &arg);			//�������� ��������� �� �������.
};

CMatrix3 operator *(FLOATTYPE scalar, const CMatrix3 &matrix);
Vec3 operator *(const Vec3 &vec, const CMatrix3 &mat);

#ifdef INLINE_BUILD
	#include "m_Matrix3.inl"
#endif

#endif
