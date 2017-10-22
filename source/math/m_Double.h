#ifndef _M_DOUBLE_H_
#define _M_DOUBLE_H_

namespace math
{
	class CDouble
	{
	public:

		FLOATTYPE mVal;

		static const FLOATTYPE epsilon;

		CDouble(void){}
		CDouble(const CDouble &arg){mVal = arg.mVal;}
		CDouble(const FLOATTYPE &arg){mVal = arg;}
		~CDouble(void){};

		CDouble &							operator = (const CDouble &arg)
		{
			if(&arg == this)
				return *this;

			mVal = arg.mVal;

			return *this;
		}

		bool								operator == (const CDouble &arg) const
		{
			return (mVal - arg.mVal) > 0 ? (mVal - arg.mVal) < epsilon : (mVal - arg.mVal) > -epsilon;
		}

		bool								operator != (const CDouble &arg) const
		{
			return !(*this == arg);
		}

		bool								operator < (const CDouble &arg) const
		{
			return mVal < arg.mVal - epsilon;
		}

		bool								operator <= (const CDouble &arg) const
		{
			return mVal < arg.mVal + epsilon;
		}

		bool								operator > (const CDouble &arg) const
		{
			return mVal > arg.mVal + epsilon;
		}

		bool								operator >= (const CDouble &arg) const
		{
			return mVal > arg.mVal - epsilon;
		}
	};
}

#endif
