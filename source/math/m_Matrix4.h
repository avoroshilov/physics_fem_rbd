#ifndef _M_MATRIX4_H_
#define _M_MATRIX4_H_

#include "all.h"

#include "m_Matrix3.h"
#include "m_Vector3.h"

class CMatrix4
{
public:
	struct
	{
		FLOATTYPE mMatrix[4][4];
	};

	CMatrix4(void);
	CMatrix4(const CMatrix4 &matrix);
	CMatrix4(const CMatrix3 &rotate, const Vec3 &translate);
	~CMatrix4(void){}

	void									MakeZero(void);
	void									MakeIdentity(void);
	void									MakeZooming(const Vec3 &origin, const FLOATTYPE &zoom);

	FLOATTYPE *								operator[](unsigned long index);
	const FLOATTYPE *							operator[](unsigned long index) const;

	CMatrix4 &								operator =(const CMatrix4 &arg);

	CMatrix4								operator *(FLOATTYPE arg) const;
	Vec3								operator *(const Vec3 &arg) const;
	CMatrix4								operator *(const CMatrix4 &arg) const;

	Vec3								RotateOnly(const Vec3 &arg) const;

	CMatrix4 &								Add(const CMatrix4 &arg);

	void									MulTranslationByIt(const Vec3 &t);//unused :)
	void									MulItByTranslation(const Vec3 &t);//unused :)
};

CMatrix4 operator *(FLOATTYPE scalar, const CMatrix4 &matrix);

#ifdef INLINE_BUILD
	#include "m_Matrix4.inl"
#endif

#endif
