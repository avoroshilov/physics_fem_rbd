#include "m_VectorN.h"
#include <math.h>

using namespace math;

CVectorN::CVectorN(void)
{
	mSize = 1;
	mVector = new FLOATTYPE[1];

	assert_t(mVector && "Error: not enough memory!");
}

CVectorN::CVectorN(unsigned long size)
{
	assert(size > 0);

	if(size == 0)
		size = 1;

	mSize = size;
	mVector = new FLOATTYPE[size];

	assert_t(mVector && "Error: not enough memory!");
}

CVectorN::CVectorN(const CVectorN &rh)
{
	mSize = 0;
	mVector = NULL;
	*this = rh;
}

CVectorN::CVectorN(const Vec3 &rh)
{
	mSize = 0;
	mVector = NULL;
	*this = rh;
}


CVectorN::~CVectorN(void)
{
	if(mVector)
		delete [] mVector;

	mVector = NULL;//For a case of FLOATTYPE-calling destructor.
}
//#pragma message("FD warning: goddamned CVectorN::Set() method wants FLOATTYPE arguments")
bool CVectorN::Set(unsigned long begin, unsigned long n, ...)
{
	assert_t(mVector);

	if(begin + n > mSize)
		return false;

	va_list argptr;
	va_start(argptr, n);
//	cout << *((FLOATTYPE* ) argptr)<<"\t";
	for(unsigned long i = begin; i < begin + n; i++)
	{
		mVector[i] = va_arg(argptr, FLOATTYPE);
//		cout << mVector[i];
	}

	va_end(argptr);

	return true;
}

FLOATTYPE CVectorN::GetLength(void) const
{
	assert_t(mVector);
	FLOATTYPE sum = 0.0;
	
	for(unsigned long i = 0; i < mSize; i++)
		sum += (mVector[i] * mVector[i]);

	return sqrt(sum);
}

FLOATTYPE CVectorN::GetSquareLength(void) const
{
	assert_t(mVector);
	
	FLOATTYPE sum = 0.0f;
	
	for(unsigned long i = 0; i < mSize; i++)
		sum += mVector[i] * mVector[i];
	
	return sum;
}

void CVectorN::MakeInverse(void)
{
	assert_t(mVector);
	
	for(unsigned long i = 0; i < mSize; i++)
		mVector[i] = -mVector[i];
	
	return;
}

FLOATTYPE CVectorN::Dot(const CVectorN &arg) const
{
	if(mSize != arg.mSize)
	{
		assert(false &&"In VectorN::Dot() sizes don't match");

		return 0.0;
	}

	FLOATTYPE sum = 0.0;
	
	for(unsigned long i = 0; i < mSize; i++)
		sum += mVector[i] * arg.mVector[i];
	
	return sum;
}



CVectorN CVectorN::ProjectOnto(const CVectorN &arg) const
{
	return arg * this->Dot(arg) / arg.GetSquareLength();
}

void CVectorN::Normalize(void)
{
	FLOATTYPE fMagnitude = this->GetLength();

	if(fMagnitude <= M_EPSILON)
		return;
	
	for(unsigned long i = 0; i< mSize; i++)
		mVector[i] /= fMagnitude;
}

CVectorN CVectorN::GetNormalized(void) const
{
	CVectorN Result(mSize);

	FLOATTYPE fMagnitude = this->GetLength();

	if(fMagnitude <= M_EPSILON)
	{
		assert(false && "Magnitude of a vector is lower than M_EPSILON!");
	//	gLog.Print("Magnitude of a vector is lower than M_EPSILON!");
		Result.MakeZero();

		return Result;
	}

	for(unsigned long i = 0; i < mSize; i++)
		Result.mVector[i] = mVector[i] / fMagnitude;

	return Result;
}

CVectorN &CVectorN::operator =(const CVectorN &arg)
{
	if(&arg == this)
		return *this;

	if(mSize != arg.GetSize())
		this->SetSize(arg.GetSize());

	for(unsigned long i = 0; i < mSize; i++)
		mVector[i] = arg.mVector[i];

	return *this;
}

#define GET_ELEMENT(a, c) ((c == 0) ? a.x : ((c == 1)? a.y : a.z ))

CVectorN &CVectorN::operator =(const Vec3 &arg)
{
	if(mSize != 3)
		this->SetSize(3);

	for(unsigned long i = 0; i < mSize; i++)
		mVector[i] = GET_ELEMENT(arg, i);

	return *this;
}

CVectorN::operator Vec3()
{
	if(mSize != 3)
	{
		assert(false && "Dimentions mismatch in Vec3 CVectorN::operator Vec3()");
		Vec3 dummy(0, 0, 0);

		return dummy;
	}

	return Vec3(mVector[0], mVector[1], mVector[2]);
}

bool CVectorN::operator ==(const CVectorN &arg) const
{
	bool IsEqual = true;
	
	for(unsigned long i = 0; i < mSize; i++)
		IsEqual = IsEqual && M_EQUAL(mVector[i], arg.mVector[i]);
	
	return IsEqual;
}

bool CVectorN::operator !=(const CVectorN &arg) const
{
	bool IsEqual = true;
	
	for(unsigned long i = 0; i < mSize; i++)
		IsEqual = IsEqual && M_EQUAL(mVector[i], arg.mVector[i]);
	
	return !IsEqual;
}


bool CVectorN::IsEqualAt(const CVectorN &arg, const FLOATTYPE &precision) const
{
	bool IsEqual = true;
	
	for(unsigned long i = 0; i < mSize; i++)
		IsEqual = IsEqual && M_EQUAL_AT(mVector[i], arg.mVector[i], precision);
	
	return IsEqual;
}

/*FLOATTYPE &CVectorN::operator [](unsigned long i)
{
//	assert_t((i < mSize) && "Index range check error.");
//	assert_t(mVector);

	return mVector[i];
}

const FLOATTYPE &CVectorN::operator [](unsigned long i) const
{
//	assert_t(i < mSize && "Index range check error.");
//	assert_t(mVector);

	return mVector[i];
}*/

CVectorN CVectorN::operator -(void) const//monadic -
{
	assert_t(mVector);
	CVectorN result(mSize);
	
	for(unsigned long i = 0;i < mSize; i++)
		result.mVector[i] = -mVector[i];
	
	return result;
}

CVectorN CVectorN::operator +(const CVectorN &arg) const
{
	assert_t((arg.mSize == mSize) && "Index range check error.");
	assert_t(mVector);
	CVectorN result(mSize);
	
	for(unsigned long i = 0; i < mSize; i++)
		result.mVector[i] = arg.mVector[i] + mVector[i];
	
	return result;
}

CVectorN CVectorN::operator -(const CVectorN &arg) const//binary -
{
	assert_t((arg.mSize == mSize) && "Index range check error.");
	assert_t(mVector);
	CVectorN result(mSize);
	
	for(unsigned long i = 0;i < mSize;i++)
		result.mVector[i] = mVector[i] - arg.mVector[i];
	
	return result;
}

CVectorN CVectorN::operator *(const FLOATTYPE &arg) const
{
	CVectorN result(mSize);

	for(unsigned long i = 0; i < mSize; i++)
		result[i] = mVector[i] * arg;

	return result;
}

CVectorN CVectorN::operator /(const FLOATTYPE &arg) const
{
	CVectorN result(mSize);

	for(unsigned long i = 0; i < mSize; i++)
		result[i] = mVector[i] / arg;

	return result;
}

CVectorN &CVectorN::operator +=(const CVectorN &arg)
{
	assert_t(arg.mSize == mSize);

	unsigned long size = mSize < arg.mSize ? mSize : arg.mSize;

	for(unsigned long i = 0; i < size; i++)
		mVector[i] += arg.mVector[i];

	return *this;
}

CVectorN &CVectorN::operator -=(const CVectorN &arg)
{
	assert_t(arg.mSize == mSize);

	unsigned long size = mSize < arg.mSize ? mSize : arg.mSize;

	for(unsigned long i = 0; i < size; i++)
		mVector[i] -= arg.mVector[i];

	return *this;
}

CVectorN &CVectorN::operator *=(const FLOATTYPE &arg)
{
	for(unsigned long i = 0; i< mSize;i++)
		mVector[i] *= arg;

	return *this;
}

CVectorN &CVectorN::operator /=(const FLOATTYPE &arg)
{
	for(unsigned long i = 0; i< mSize;i++)
		mVector[i] /= arg;

	return *this;
}
/*CVectorN operator *(const float &scalar, const CVectorN &vector)
{
	return vector * scalar;
}*/
void CVectorN::MakeZero(void)
{
	for(unsigned long i = 0; i < mSize; i++)
		mVector[i] = 0.0;
}

bool CVectorN::SetSize(unsigned long size)
{
	assert(size > 0);

	if(size == 0)
		size = 1;

	FLOATTYPE *OldVector = mVector;

	mVector = new FLOATTYPE[size];

	assert_t(mVector && "Error: not enough memory!");

	if(!mVector)
	{
		mVector = OldVector;

		return false;
	}

	//Copy old contents.

	unsigned long minsize = size < mSize ? size : mSize;

	for(unsigned long i = 0; i < minsize; i++)
		mVector[i] = OldVector[i];

	delete [] OldVector;

	mSize = size;

	return true;
}

unsigned long CVectorN::GetSize(void) const
{
	return mSize;
}

void CVectorN::SwapElements(unsigned long i, unsigned long j)
{
	assert_t(i < mSize && "Array range check error (i).");
	assert_t(j < mSize && "Array range check error (j).");

	if(i == j)
		return;

	FLOATTYPE t    = mVector[i];
	mVector[i] = mVector[j];
	mVector[j] = t;
}

void CVectorN::Log(void)
{
	char str[1024];
	str[0] = '\0';

	for(unsigned long r = 0; r < mSize; r++)
	{
		char s[20];
		sprintf(s, "%.4f ", mVector[r]);
		strcat(str, s);
	}

	gLog.Print("%s", str);
}
