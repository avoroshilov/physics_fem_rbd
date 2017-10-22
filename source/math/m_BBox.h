#ifndef _M_BBOX_H_
#define _M_BBOX_H_

#include "m_Vector3.h"

namespace math
{
	class CAABBox
	{
	public:

		Vec3							mPos;
		Vec3							mDims;//half width, half height, half depth.

		CAABBox(void)
		{
			mPos.Set(0, 0, 0);
			mDims.Set(1, 1, 1);
		}

		CAABBox(const CAABBox &bbox)
		{
			if(bbox.mDims.x <= 0)
				message("die");

			assert_t(bbox.mDims.x > 0);
			assert_t(bbox.mDims.y > 0);
			assert_t(bbox.mDims.z > 0);

			mPos  = bbox.mPos;
			mDims = bbox.mDims;
		}

		CAABBox(Vec3 pos, Vec3 dims)
		{
			assert_t(dims.x > 0);
			assert_t(dims.y > 0);
			assert_t(dims.z > 0);

			mPos  = pos;
			mDims = dims;
		}

		~CAABBox(void){}

		void								Set(Vec3 pos, Vec3 dims)
		{
			assert_t(dims.x > 0);
			assert_t(dims.y > 0);
			assert_t(dims.z > 0);

			mPos = pos;
			mDims = dims;
		}

		bool								Contains(const CAABBox &ref, FLOATTYPE edging = 0.0)
		{
			assert_t(edging < mDims.x);
			assert_t(edging < mDims.y);
			assert_t(edging < mDims.z);

			//test x
			if(ref.mPos.x - ref.mDims.x <= mPos.x - mDims.x + edging)
				return false;

			if(ref.mPos.x + ref.mDims.x >= mPos.x + mDims.x - edging)
				return false;

			//test y
			if(ref.mPos.y - ref.mDims.y <= mPos.y - mDims.y + edging)
				return false;

			if(ref.mPos.y + ref.mDims.y >= mPos.y + mDims.y - edging)
				return false;

			//test z
			if(ref.mPos.z - ref.mDims.z <= mPos.z - mDims.z + edging)
				return false;

			if(ref.mPos.z + ref.mDims.z >= mPos.z + mDims.z - edging)
				return false;

			return true;
		}
	};

	class CRangedAABBox
	{
	public:

		Vec3							mFrom;
		Vec3							mTo;

		CRangedAABBox(void)
		{
			mFrom.Set(0, 0, 0);
			mTo.Set(0, 0, 0);
		}

		CRangedAABBox(const CRangedAABBox &arg)
		{
			mFrom = arg.mFrom;
			mTo = arg.mTo;
		}

		~CRangedAABBox(void){}

		bool								Contains(const CAABBox &ref, FLOATTYPE edging = 0.0)
		{
			//test x
			if(ref.mPos.x - ref.mDims.x <= mFrom.x + edging)
				return false;

			if(ref.mPos.x + ref.mDims.x >= mTo.x - edging)
				return false;

			//test y
			if(ref.mPos.y - ref.mDims.y <= mFrom.y + edging)
				return false;

			if(ref.mPos.y + ref.mDims.y >= mTo.y - edging)
				return false;

			//test z
			if(ref.mPos.z - ref.mDims.z <= mFrom.z + edging)
				return false;

			if(ref.mPos.z + ref.mDims.z >= mTo.z - edging)
				return false;

			return true;
		}
	};
}

#endif
