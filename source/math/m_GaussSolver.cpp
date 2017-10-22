#include "m_GaussSolver.h"

#include "m_Const.h"

#include "helpers/debug.h"

#pragma message(__FILE__" : warning: Turn on assert macros!")

//

using namespace math;

bool CGaussSolver::SetSize(unsigned long size)
{
	return mMatrix.SetSize(size) && mVector.SetSize(size) && mResult.SetSize(size);
}

bool CGaussSolver::Solve(void)
{
	assert( (mMatrix.GetSize() == mVector.GetSize()) && (mMatrix.GetSize() == mResult.GetSize()) );

	if( (mMatrix.GetSize() != mVector.GetSize()) || (mMatrix.GetSize() != mResult.GetSize()) )
		return false;

//	mResult.MakeZero();
	const unsigned long size = mMatrix.GetSize();
	//вектор индексов. Поскольку в этом алгоритме могут меняться местами СТОБЦЫ,
	//indices[номер столбца] = номер переменной, стоящей над этим столбцом
	unsigned long* indices = new unsigned long[size];
	
	for(unsigned long i = 0; i<size;i++)
		indices[i] = i;

//	com::CArray<bool> OmittedVar;
//	OmittedVar.SetUsedSize(size);
//	com::CArray<unsigned long> StartVar;
//	StartVar.SetUsedSize(size);
//	unsigned long shift = 0;
/*	for(unsigned long i = 0; i< size; i++)
	{
		StartVar[i]= i;
		OmittedVar[i]= false;
		mResult[i]= 0.0f;
	}*/

	//-----------------------------------------------------------------------------
	// Первый проход
	//-----------------------------------------------------------------------------
	//внешний цикл - для каждой переменной. В этом обязательно разрешается
	//1 переменная
	unsigned long var = 0;

	for(var = 0;var<size;var++)
	{
		//В среднем цикле идёт поиск столбца (из оставшихся) с хотя бы одним
		//ненулевым элементом
		unsigned long col;
		for(col = var;col<size;col++)
		{
			//Во внутренним цикле идёт поиск строки с ненулевым элементом в даннйо строке

			unsigned long row;
			for(row = var;row<size;row++)
			{
				//если такой элемент найден...
				if(mMatrix(row, col)!= 0)
				{
					//и если он не первый возможный (т.е. стоящей в строке,
					//соответствующей номеру разрешаемой переменной), то надо
					//поменять сстроки местами
					if(row!= var)
					{
						FLOATTYPE temp;

						for(unsigned long i = var;i<size;i++)
						{
							temp = mMatrix(var, i);
							mMatrix(var, i)= mMatrix(row, i);
							mMatrix(row, i) = temp;
						}

						temp = mVector[var];
						mVector[var]= mVector[row];
						mVector[row]= temp;
					}
				//	Если найден ненулевой элемент (и строка с ним переставлена
				// в нужное место, то дальше просматривать элементы других строк
					//данного столбца ненужно

					break;
				}
			}
			//Выполненение этого условия означает, что в данном столбце был
			//найден ненулевой элемент
			if(row!= size)
			{
				//если этот столбец - не первый возможный, т.е. были полностью нулевые
				//столбцы, нужно это тстолбец поставить на нужное место, чтобы
				//результатирующая матрица была всё равно треугольной
				if(col!= var)
				{
					FLOATTYPE temp;
					for(unsigned long i = 0; i<size;i++)
					{
						temp = mMatrix(i, var);
						mMatrix(i, var) = mMatrix(i, col);
						mMatrix(i, col) = temp;
					}
					//За счёт этого массива устанавливаем связь между фактическим столбцом
					//и преременной, за которую он отвечает.
					unsigned long ltemp = indices[var];
					indices[var] = indices[col];
					indices[col]= ltemp;
				}
				//если найден столбец с ненулевым элементом, дальше не ищем
				break;
			}
		}
		//это условие означает, что столбец с ненулевым элементом НЕ найден,
		//то есть мы имеем треугольную матрицу, но все стркоки включая и ниже текущей
		//перменной НУЛЕВЫЕ, т.е. мы имеем либо множество, либо ноль решений.
		//в таком случае уже текущую прееменную нельязя решить. Прерываем внешний цикл
		if(col == size)
		{
			break;
		}
		//Если мы здесь, то по адрессу (var, var) стоит ненулевой элемент
		//Разрешаем для этой перменной матрицу.
		//FLOATTYPE div = 1/mMatrix(var, var);
		//шаг 1. Деление всей строк ина разрешающий элемент
		FLOATTYPE d = mMatrix(var, var);
		mVector[var] = mVector[var] / d;

		for(unsigned long i = var; i<size;i++)
			mMatrix(var, i) = mMatrix(var, i) / d;

		//шаг 2. Вычитаение данной строки ( с коеффициентом) из всех строк ниже данной.
		for(unsigned long i = var+1; i<size;i++)
			mVector[i] -= mVector[var] * mMatrix(i, var);

		for(unsigned long j = var+1; j<size;j++)
		{
			FLOATTYPE sub = mMatrix(j, var);
			for(unsigned long i = var; i<size;i++)
				mMatrix(j, i) -= mMatrix(var, i)*sub;
		}
	}
	//*******************************
	//Проверка наличия решения
	//*******************************

	//это условие значит, что число разрешённых элементов меньше изначального
	//числа переменных. Т.е. внизу матрицы есть нулевые строки. Если элементы
	//вектора, соотвтетсвующие этим строкам ненулевые, то решений нет, иначе
	//решений множество.
	if(var != size)
	{
		for(unsigned long i = var; i<size; i++)
		{
			//Используем мягко стравнение, потмоу что "ошибка в эту сторону"
			//не так страшна, как отсутствие решения вообще
			if(!M_EQUAL_AT(mVector[i], 0, 0.000001f))
			{
				delete[] indices;

				return false;
			}
		}
	}

	//******************************
	//Обратный ход гаусса
	//******************************

	//здесь для начала Вар - первая переменная, для которой нельяз найти решение
	//принимаем всё свободные переменные за нули.
	for(unsigned long i = var;i>0;i--)
	{
	//	v[i] /= m[i*dim + i];
		//Прочистка столбца над данной перменной
		//-1 везде используется, потому что переменные цикла больше 1 своего
		//"физического смысла", для того, чтоб цикл никогда не уходил в минус -
		//таки беззнаковость кругом.
		for(unsigned long j = i-1; j>0;j--)
		{

			mVector[j-1] -= mMatrix(j-1, i-1)*mVector[i-1];
		}
	}
	//**************************************
	//Восстанвление порядка переменных
	//**************************************
	for(unsigned long i = 0; i<size;i++)
	{
		mResult[indices[i]]= mVector[i];
	}

	delete[] indices;

	return true;
}

bool CGaussSolver::Solve0(void)
{
	assert( (mMatrix.GetSize() == mVector.GetSize()) && (mMatrix.GetSize() == mResult.GetSize()) );

	if( (mMatrix.GetSize() != mVector.GetSize()) || (mMatrix.GetSize() != mResult.GetSize()) )
		return false;
//	mResult.MakeZero();
	const unsigned long size = mMatrix.GetSize();
	//вектор индексов. Поскольку в этом алгоритме могут меняться местами СТОБЦЫ,
	//indices[номер столбца] = номер переменной, стоящей над этим столбцом
	unsigned long* indices = new unsigned long[size];
	
	for(unsigned long i = 0; i<size;i++)
		indices[i] = i;

//	com::CArray<bool> OmittedVar;
//	OmittedVar.SetUsedSize(size);
//	com::CArray<unsigned long> StartVar;
//	StartVar.SetUsedSize(size);
//	unsigned long shift = 0;
/*	for(unsigned long i = 0; i< size; i++)
	{
		StartVar[i]= i;
		OmittedVar[i]= false;
		mResult[i]= 0.0f;
	}*/

	//-----------------------------------------------------------------------------
	// Первый проход
	//-----------------------------------------------------------------------------
	//внешний цикл - для каждой переменной. В этом обязательно разрешается
	//1 переменная
	unsigned long var = 0;
	for(var = 0;var<size;var++)
	{
		//В среднем цикле идёт поиск столбца (из оставшихся) с хотя бы одним
		//ненулевым элементом
		unsigned long col;
		for(col = var;col<size;col++)
		{
			//Во внутренним цикле идёт поиск строки с ненулевым элементом в даннйо строке

			unsigned long row;
			for(row = var;row<size;row++)
			{
				//если такой элемент найден...
				if(mMatrix(row, col)!= 0)
				{
					//и если он не первый возможный (т.е. стоящей в строке,
					//соответствующей номеру разрешаемой переменной), то надо
					//поменять сстроки местами
					if(row!= var)
					{
						FLOATTYPE temp;
						for(unsigned long i = var;i<size;i++)
						{
							temp = mMatrix(var, i);
							mMatrix(var, i)= mMatrix(row, i);
							mMatrix(row, i) = temp;
						}
						temp = mVector[var];
						mVector[var]= mVector[row];
						mVector[row]= temp;
					}
				//	Если найден ненулевой элемент (и строка с ним переставлена
				// в нужное место, то дальше просматривать элементы других строк
					//данного столбца ненужно

					break;
				}
			}
			//Выполненение этого условия означает, что в данном столбце был
			//найден ненулевой элемент
			if(row!= size)
			{
				//если этот столбец - не первый возможный, т.е. были полностью нулевые
				//столбцы, нужно это тстолбец поставить на нужное место, чтобы
				//результатирующая матрица была всё равно треугольной
				if(col!= var)
				{
					FLOATTYPE temp;
					for(unsigned long i = 0; i<size;i++)
					{
						temp = mMatrix(i, var);
						mMatrix(i, var) = mMatrix(i, col);
						mMatrix(i, col) = temp;
					}
					//За счёт этого массива устанавливаем связь между фактическим столбцом
					//и преременной, за которую он отвечает.
					unsigned long ltemp = indices[var];
					indices[var] = indices[col];
					indices[col]= ltemp;
				}
				//если найден столбец с ненулевым элементом, дальше не ищем
				break;
			}
		}
		//это условие означает, что столбец с ненулевым элементом НЕ найден,
		//то есть мы имеем треугольную матрицу, но все стркоки включая и ниже текущей
		//перменной НУЛЕВЫЕ, т.е. мы имеем либо множество, либо ноль решений.
		//в таком случае уже текущую прееменную нельязя решить. Прерываем внешний цикл
		if(col == size)
		{
			break;
		}
		//Если мы здесь, то по адрессу (var, var) стоит ненулевой элемент
		//Разрешаем для этой перменной матрицу.
		FLOATTYPE div = 1/mMatrix(var, var);
		//шаг 1. Деление всей строк ина разрешающий элемент
		mVector[var] *= div;

		for(unsigned long i = var; i<size;i++)
			mMatrix(var, i) *= div;

		//шаг 2. Вычитаение данной строки ( с коеффициентом) из всех строк ниже данной.
		for(unsigned long i = var+1; i<size;i++)
			mVector[i] -= mVector[var] * mMatrix(i, var);

		for(unsigned long j = var+1; j<size;j++)
		{
			FLOATTYPE sub = mMatrix(j, var);
			for(unsigned long i = var; i<size;i++)
				mMatrix(j, i) -= mMatrix(var, i)*sub;
		}
	}
	//*******************************
	//Проверка наличия решения
	//*******************************

	//это условие значит, что число разрешённых элементов меньше изначального
	//числа переменных. Т.е. внизу матрицы есть нулевые строки. Если элементы
	//вектора, соотвтетсвующие этим строкам ненулевые, то решений нет, иначе
	//решений множество.
	if(var != size)
	{
		for(unsigned long i = var; i<size; i++)
		{
			//Используем мягко стравнение, потмоу что "ошибка в эту сторону"
			//не так страшна, как отсутствие решения вообще
			if(!M_EQUAL_AT(mVector[i], 0, 0.000001f))
			{
				delete[] indices;

				return false;
			}
		}
	}

	//******************************
	//Обратный ход гаусса
	//******************************

	//здесь для начала Вар - первая переменная, для которой нельяз найти решение
	//принимаем всё свободные переменные за нули.
	for(unsigned long i = var;i>0;i--)
	{
	//	v[i] /= m[i*dim + i];
		//Прочистка столбца над данной перменной
		//-1 везде используется, потому что переменные цикла больше 1 своего
		//"физического смысла", для того, чтоб цикл никогда не уходил в минус -
		//таки беззнаковость кругом.
		for(unsigned long j = i-1; j>0;j--)
		{

			mVector[j-1] -= mMatrix(j-1, i-1)*mVector[i-1];
		}
	}
	//**************************************
	//Восстанвление порядка переменных
	//**************************************
	for(unsigned long i = 0; i<size;i++)
	{
		mResult[indices[i]]= mVector[i];
	}

	delete[] indices;

	return true;
}

bool CGaussSolver::Solve1(void)
{
	assert( (mMatrix.GetSize() == mVector.GetSize()) && (mMatrix.GetSize() == mResult.GetSize()) );

	if( (mMatrix.GetSize() != mVector.GetSize()) || (mMatrix.GetSize() != mResult.GetSize()) )
		return false;

	const unsigned long size = mMatrix.GetSize();

	//-----------------------------------------------------------------------------
	// Prepare the matrix
	//-----------------------------------------------------------------------------

	for(unsigned long n = 0; n < size - 1; n++)
	{
		//-----------------------------------------------------------------------------
		// Find a base row
		//-----------------------------------------------------------------------------

		for(unsigned long i = n; i < size; i++)
		{
			//TODO: replace here with abs.
			if(mMatrix(i, n) < -0.000001f || mMatrix(i, n) > 0.000001f)
			{
				mMatrix.SwapRows(n, i);
				mVector.SwapElements(n, i);

				break;
			}

			if(i == size - 1)//Bad matrix.
				return false;
		}

		//-----------------------------------------------------------------------------
		// A pass
		//-----------------------------------------------------------------------------

		for(unsigned long i = n + 1; i < size; i++)
		{
			if(mMatrix(n, n) > -0.000001f && mMatrix(n, n) < 0.000001f)
				return false;

			FLOATTYPE k = mMatrix(i, n) / mMatrix(n, n);

			for(unsigned long j = n + 1; j < size; j++)
				mMatrix(i, j) -= mMatrix(n, j) * k;

			mVector[i] -= mVector[n] * k;
		}
	}

	//-----------------------------------------------------------------------------
	// Extract results
	//-----------------------------------------------------------------------------

	for(unsigned long i = size - 1; ; i--)
	{
		FLOATTYPE res;

		if(mMatrix(i, i) > -0.0001f && mMatrix(i, i) < 0.0001f)
		{
			mResult[i] = 1.0f;
		}
		else
		{
			res = mVector[i];

			for(unsigned long j = i + 1; j < size; j++)
				res -= mMatrix(i, j) * mResult[j];

			mResult[i] = res / mMatrix(i, i);
		}

		if(i == 0)
			break;
	}

	return true;
}