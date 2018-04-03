#pragma once
#include "Setups.h"
#include <initializer_list>
#include <vector>
//using namespace std;
#define EQUAL(X,Y) ( ABS((X) - (Y)) <= 1e-6 ? (Y) : (X) )
#define ABS(X) ( (X) < 0 ? -(X) : (X) )

template <typename T>

class Matrix
{ 
protected:
	
	

public:
	std::vector<T> data;
	size_t m, n;

	void swap_rows(const size_t &First, const size_t &Second);
	void simple_converting();
	//����� �������
	void Print();
	//***********************
	//������������

	Matrix();
	/*Matrix(Matrix<T> _data);*/
	Matrix(Matrix<T>& _data);
	Matrix(const std::vector<std::vector<T>> &matr);
	Matrix(T *A,const size_t _m, const size_t _n);
	Matrix(const size_t &NumberOfRows,const size_t &NumberOfColumns);
	//Matrix(bool IdentityMatrix, size_t Size);//Size - ������ ���������� �������, Indentity - �������� �� ��������� ��������
	static Matrix<T> MakeIdentity(size_t Size);

	//***********************
	//��������� ������

	T& operator () (const size_t &_M, const size_t &_N);
	const T& operator () (const size_t &_M,const size_t &_N) const;
	Matrix<T> operator = (const Matrix<T> &_data);
	Matrix<T> operator + (const Matrix<T> &_Var);
	Matrix<T> operator * (const T &_Value);			//��������� �� �����
	Matrix<T> operator * (const Matrix<T> &_Var);	//��������� ���������
	Matrix<T> operator *= (const Matrix<T> &_Var);	//������������ ���������
	Matrix<T> operator - (const Matrix<T> &_Var);
	//***********************
	Matrix<T> transpose();
	float determinant();
	Matrix<T> invers();
	Matrix<T> augment(const Matrix<T> &B);
	//***********************
};
template<typename T>
inline void Matrix<T>::swap_rows(const size_t & First,const size_t & Second)
{
	for (size_t j = 0; j < n; j++)
	{
		T temp = (*this)(Second, j);
		(*this)(Second, j) = (*this)(First, j);
		(*this)(First, j) = temp;
	}
}
template<typename T>
inline void Matrix<T>::simple_converting()
{
	for (size_t i = 1; i < m; i++)
	{
		T top = (*this)(i, 0);
		T down = (*this)(0, 0);
		for (size_t j = 0; j < n; j++)
		{
			auto ANS = (*this)(i, j) - (top / down)*((*this)(0,j));
			(*this)(i, j) = ANS;
		}
	}
	for (size_t i = 0; i < m; i++)
	{
		this->data.erase(data.begin());
	}
	m--;
	n--;
	for (size_t i = 0; i < m; i++)
	{
		this->data.erase(data.begin() + i*m);
	}
}
//����� �������. ����� �������
template<typename T>
inline void Matrix<T>::Print()
{
	for (size_t i = 0; i < m; i++)
	{
		cout << "\n ";
		for (size_t j = 0; j < n; j++)
			cout << data[i*n + j] << " ";
	}
	cout << endl;
}
//****************************

template<typename T>
inline Matrix<T>::Matrix()
	:m(4), n(4)
{
	data.reserve(m*n);
	for (size_t i = 0; i < m; i++)
	{vector<T>
		for (size_t j = 0; j < n; j++)
		{
			data.emplace_back(0.0);
		}
	}
}

template<typename T>
inline Matrix<T>::Matrix(Matrix<T>&  _data)
	:m(_data.m),n(_data.n),data(_data.data)
{
}

template<typename T>
inline Matrix<T>::Matrix(const std::vector<std::vector<T>>& matr):m(matr.size()),n(matr[0].size())
{
	data.reserve(m*n);
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			data.emplace_back(matr[i].at(j));
		}
	}

}

template<typename T>
inline Matrix<T>::Matrix(T *A, const size_t _m, const size_t _n)
	: m(_m), n(_n)
{
	for (size_t i = 0; i < m; i++)
	{
		for (size_t i = 0; i < n; i++)
		{
			this->(i, j) = A[i][j];
		}
	}
}

template<typename T>
inline Matrix<T>::Matrix(const size_t &NumberOfRows, const size_t &NumberOfColumns)
	: m(NumberOfRows), n(NumberOfColumns)
{
	data.reserve(m*n);
	for (size_t i = 0; i < m; i++)
	{
		for (size_t j = 0; j < n; j++)
		{
			data.emplace_back(0.0);
		}
	}
}

template<typename T>
inline Matrix<T> Matrix<T>::MakeIdentity(size_t Size)
{
	auto m = Matrix<T>(Size, Size);
	for (int i = 0; i < Size; i++)
	{
		m(i, i) = 1.0;
	}
	return m;
}


template<typename T>
inline Matrix<T> Matrix<T>::operator=(const Matrix<T> &_data)
{
	m = _data.m;
	n = _data.n;
	data = _data.data;
	return *this;
}

template<typename T>
inline Matrix<T> Matrix<T>::operator+(const Matrix<T> &_Var)
{
	Matrix<T> Result;
	Result.m = _Var.m;
	Result.n = _Var.n;
	Result.data.reserve(m*n);
	if ((_Var.m != m) || (_Var.n != n))	return Matrix<T>();
	else
	{
		size_t index;
		for (auto Iter = _Var.data.begin(); Iter != _Var.data.end(); Iter++)
		{
			index = Iter - _Var.data.begin();
			Result.data.at(index) = data.at(index) + _Var.data.at(index);
		}
		return Result;
	}
}

template<typename T>
inline Matrix<T> Matrix<T>::operator*(const T &_Value)
{
	Matrix<T> Result(m,n);
	size_t index;
	for (auto Iter = data.begin(); Iter != data.end(); Iter++)
	{
		index = Iter - data.begin();
		Result.data.at(index) = data.at(index) * _Value;
	}
	return Result;
}

template<typename T>
inline Matrix<T> Matrix<T>::operator*(const Matrix<T> &_Var)
{
	if (n != _Var.m)	return Matrix<T>();
	else
	{
		Matrix<T> Result(m,_Var.n);
		for (size_t i = 0; i < Result.m; i++)
		{
			for (size_t j = 0; j < Result.n; j++)
			{
				for (size_t k = 0; k < n; k++)
				{
					Result(i, j) += (*this)(i, k)*_Var(k, j);
					Result(i, j) = EQUAL(Result(i, j), 1.0);
					Result(i, j) = EQUAL(Result(i, j), 0.0);
				}
			}
		}
		return Result;
	}
}

template<typename T>
inline Matrix<T> Matrix<T>::operator *= (const Matrix<T> &_Var)
{
	if ((m != _Var.m)||(n != _Var.n)) return Matrix<T>();
	else
	{
		Matrix<T> Result(m, n);
		for (size_t i = 0; i < m; i++)
		{
			for (size_t j = 0; j < n; j++) Result(i, j) = (*this)(i, j)*_Var(i, j);
		}
		return Result;
	}
	
}

template<typename T>
inline Matrix<T> Matrix<T>::operator-(const Matrix<T>& _Var)
{
	Matrix<T> result(_Var.m, _Var.n);

	for (size_t i = 0; i < result.m; i++)
	{
		for (size_t j = 0; j < result.n; j++)
		{
			result(i, j) = (*this)(i,j) - _Var(i, j);
		}
	}
	return result;
}

template<typename T>
inline Matrix<T> Matrix<T>::transpose()
{
	Matrix<T> Result(n,m);
	for (size_t i = 0; i < m; i++)
	{
		for (size_t j = 0; j < n; j++)
		{
			Result(j, i) = (*this)(i, j);
		}
	}
	return Result;
}

template<typename T>
inline float Matrix<T>::determinant()
{
	if (m != n) return 1;
	else
	{
		int change_sign = 1;
		Matrix<T> temp = (*this);
		if (temp(0, 0) == 0)
		{
			size_t i = 0;
			do
			{
				i++;
				if (i == temp.m)	return 0.0;
			} while (temp(i,0)==0);
			size_t first_row = 0;
			temp.swap_rows(first_row, i);	//������ ������� ������ ������ � i-��, ����� ������ ������� ��� ���������.
			change_sign *= -1;
		}
		
		double result = 1;	//������ ������� �������, ������� ��������� �� �������
		for (size_t i = 0; i < m - 2; i++)
		{
			result *= temp(0, 0);
			temp.simple_converting();
			if (temp(0, 0) == 0)
			{
				size_t i = 0;
				do
				{
					i++;
					if (i == temp.m)	return 0.0;
				} while (temp(i, 0) == 0);
				size_t first_row = 0;
				temp.swap_rows(first_row, i);	//������ ������� ������ ������ � i-��, ����� ������ ������� ��� ���������.
				change_sign *= -1;
			}
		}
		result *= (temp(0, 0)*temp(1, 1) - temp(1, 0)*temp(0, 1))*change_sign;
		return result;
	}
}

template<typename T>
inline Matrix<T> Matrix<T>::invers()
{
	if (m != n) return Matrix<T>();
	else
	{
		Matrix<T> temp((*this));
		Matrix<T> result = MakeIdentity(temp.m);
		if ((*this)(0,0) == 0)
		{
			size_t number = 0;
			do
			{
				number++;
				if (number == temp.m) return Matrix<T>();	//����������, ��� ���������� � ���� ������
			} while (temp(number, 0) == 0);
			size_t first_row = 0;
			temp.swap_rows(first_row, number);
			result.swap_rows(first_row, number);
		}
		for (size_t i = 0; i < temp.m; i++)
		{
			T vip_element = temp(i, i);
			for (size_t j = 0; j < temp.m; j++)
			{
				temp(i, j) = EQUAL(temp(i, j) / vip_element,0.0);
				result(i, j) = EQUAL(result(i, j) / vip_element,0.0);
			}
			for (size_t row = 0; row < temp.m; row++)
			{
				if (row != i)
				{
					if (temp(row, i) != 0)
					{
						T first_element = temp(row, i);
							for (size_t k = 0; k < temp.m; k++)
							{
								temp(row, k) = EQUAL(temp(row, k) - temp(i, k)*first_element,0.0);
								result(row, k) = EQUAL(result(row, k) - result(i, k)*first_element,0.0);
							}
					}
					else continue;
				}

			}
		}
		return result;
	}
}

template<typename T>
inline Matrix<T> Matrix<T>::augment(const Matrix<T>& B)
{
	if ((m != B.m)) throw(exception());
	Matrix<T> A(m, B.n + n);
	A.data.reserve(m*n + B.m * B.n);
	for (size_t i = 0; i < m; i++)
	{
		for (size_t j = 0; j < n + B.n; j++)
		{
			if (j >= n) A(i,j) = B(i, j - n);
			else	A(i,j) = (*this)(i, j);
		}
	}
	return A;
}

template<typename T>
inline T& Matrix<T>::operator () (const size_t &_M,const size_t &_N)
{
	return data[_M*n + _N];
}

template<typename T>
inline const T& Matrix<T>::operator () (const size_t &_M,const size_t &_N) const
{
	return data[_M*n + _N];
}