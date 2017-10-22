#ifndef _CAMERA_H_
#define _CAMERA_H_

#include "all.h"
#include "math/m_Vector3.h"

const Vec3 Up(0.0, 1.0, 0.0);

class CCamera 
{
public:
	Vec3 mPos, mViewDir, mStrafeVec;

	bool mbReturnMmPos;
	long mOldX, mOldY;
	long mMiddleX, mMiddleY;
	
	float mSens;

	CCamera();
	
	void SetSensivity(float sens);
	float GetSensivity() const;

	void SetReturnMouseToCenter(bool bm);
	bool GetReturnMouseToCenter() const;

	void SetScreenCenterCoords(long middlex, long middley);
	void GetSavedMouseCoords(long &x, long &y) const;
	void Setup(const Vec3 &pos, const Vec3 &target);
	
	const Vec3 &GetPos() const;
	const Vec3 &GetViewdir() const;
	const Vec3 &GetStrafedir() const;

	void OnMouse(long x, long y);

	void Strafe(float speed);
	void MoveXZ(float speed);
	void MoveForward(float speed);
	void MoveUp(float speed);
};

#endif