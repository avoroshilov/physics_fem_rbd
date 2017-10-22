#include "camera.h"

#include "math/m_Quaternion.h"

CCamera::CCamera()
{
	mOldX = -1;
	mOldY = -1;
	
	mMiddleX = -1;
	mMiddleY = -1;

	mbReturnMmPos = true;

	mSens = 1.0f;
}
void CCamera::SetScreenCenterCoords(long middlex, long middley)
{
	mMiddleX = middlex;
	mMiddleY = middley;
}
		
void CCamera::SetSensivity(float sens)
{
	mSens = sens;
}

float CCamera::GetSensivity() const
{
	return mSens;
}

void CCamera::SetReturnMouseToCenter(bool bm)
{
	mbReturnMmPos = bm;
}

bool CCamera::GetReturnMouseToCenter() const
{
	return mbReturnMmPos;
}

void CCamera::GetSavedMouseCoords(long &x, long &y) const
{
	x = mOldX;
	y = mOldY;
}

void CCamera::Setup(const Vec3 &pos, const Vec3 &target)
{
	mPos = pos;
	mViewDir = (target - mPos).GetNormalized();
	mStrafeVec = mViewDir.Cross(Up);
	mStrafeVec.Normalize();
}

const Vec3 &CCamera::GetPos() const
{
	return mPos;
}

const Vec3 &CCamera::GetViewdir() const
{
	return mViewDir;
}

const Vec3 &CCamera::GetStrafedir() const
{
	return mStrafeVec;
}

void CCamera::OnMouse(long x, long y)
{
	if(mOldX == -1 || mOldY == -1)
	{
		mOldX = x;
		mOldY = y;

		return;
	}

	if ((x == mOldX) && (y == mOldY))
		return;

	float angleY = ((mOldX - x)) / 500.0f * mSens;		
	float angleZ = ((mOldY - y)) / 500.0f * mSens;		

	if (mbReturnMmPos)
	{
		mOldX = mMiddleX;
		mOldY = mMiddleY;
	}
	else
	{
		mOldX = x;
		mOldY = y;
	}

	Vec3 vAxis = mViewDir.Cross(Up);
	vAxis.Normalize();
	CQuaternion q(vAxis, angleZ);
	Vec3 tmpView = q * mViewDir;
	tmpView.Normalize();

	Vec3 tmpVec = tmpView.Cross(Up);

	if(tmpVec.Dot(vAxis) > 0.0f)
	{
		mViewDir = tmpView;
	}

	CQuaternion quat(Up, angleY);
	mViewDir = quat * mViewDir;

	mStrafeVec = mViewDir.Cross(Up);
	mStrafeVec.Normalize();
}

void CCamera::Strafe(float speed)
{
	mPos.x += mStrafeVec.x * speed;
	mPos.z += mStrafeVec.z * speed;
}

void CCamera::MoveXZ(float speed)
{
	Vec3 vTemp(mViewDir.x, 0.0f, mViewDir.z);
	
	if(vTemp.GetLength() < M_EPSILON)
	{
		vTemp.Set(1.0f, 0.0f, 0.0f);
	}
	else
	{
		vTemp.Normalize();
	}

	mPos += vTemp * speed;
}

void CCamera::MoveForward(float speed)
{
	mPos += mViewDir * speed;
}

void CCamera::MoveUp(float speed)
{
	mPos.y += speed;
}

