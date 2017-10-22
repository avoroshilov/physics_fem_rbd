#include "m_SetOfIndices.h"

#include "helpers/debug.h"

using namespace math;

CSetOfIndices::CSetOfIndices(unsigned long MaxSize)
{
	mSize = 0;
	mMaxSize = MaxSize;
	mIndices = new unsigned long[MaxSize];

	assert_t(mIndices);
}

CSetOfIndices::~CSetOfIndices(void)
{
	if(mIndices)
	{
		delete [] mIndices;
		mIndices = NULL;
	}
}

unsigned long &CSetOfIndices::operator [](unsigned long i)
{
	assert_t(i < mSize);

	return mIndices[i];
}

const unsigned long &CSetOfIndices::operator [](unsigned long i) const
{
	assert_t(i < mSize);

	return mIndices[i];
}

void CSetOfIndices::Add(unsigned long index)
{
	assert(mSize <= mMaxSize);
	assert(index < mMaxSize);

	if(mSize >= mMaxSize)//The list is full.
		return;

	if(index >= mMaxSize)//Index is out of acceptable range.
		return;

	for(unsigned long i = 0; i < mSize; i++)
	{
		if(mIndices[i] == index)
			return;
	}

	mIndices[mSize] = index;
	mSize++;
}

void CSetOfIndices::Delete(unsigned long index)
{
	for(unsigned long i = 0; i < mSize; i++)
	{
		if(mIndices[i] == index)
		{
			mIndices[i] = mIndices[mSize - 1];
			mSize--;

			break;
		}
	}
}

void CSetOfIndices::Empty(void)
{
	mSize = 0;
}

bool CSetOfIndices::Contains(unsigned long index)
{
	for(unsigned long i = 0; i < mSize; i++)
	{
		if(mIndices[i] == index)
			return true;
	}

	return false;
}
