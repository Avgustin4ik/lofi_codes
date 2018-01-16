#pragma once

template <typename T, size_t M, size_t N>
void swap_rows(T (&A)[M][N], const size_t &first, const size_t &second)
{
	T temp_raw[N];
	for (size_t j = 0; j < N; j++)
	{
		temp_raw[j] = A[first][j];
	}
	for (size_t j = 0; j < N; j++)
	{
		A[first][j] = A[second][j];
	}
	for (size_t j = 0; j < N; j++)
	{
		A[second][j] = temp_raw[j];
	}
}
template <typename T, size_t M, size_t N>
size_t find_max_of_selected_column(T(&A)[M][N], const size_t &_j, const size_t &i_top)
{
	T max;
	size_t i_max;
	max = fabs(A[i_top][_j]);
	i_max = i_top;
	for (size_t i = i_top; i < M; i++)
	{
		if (fabs(A[i][_j]) > max)
		{
			max = A[i][_j];
			i_max = i;
		}
	}
	return i_max;
}
template <typename T, size_t M, size_t N>
T advanced_summ_of_row(T(&A)[M][N], T(&X)[M][1], const size_t &_i)
{
	T summ = 0;
	for (size_t j = _i + 1; j < N; j++)
	{
		summ = summ + A[_i][j] * X[j][0];
	}
	return summ;
}
template <typename T, size_t M, size_t N>
void method_Gauss_SLAU(T(&A)[M][N], T(&B)[M][1], T(&X)[M][1])
{
	//*********direct order********
	for (size_t k = 0; k < M-1; k++)
	{
		size_t index_max = find_max_of_selected_column(A, k, k);
		swap_rows(A, index_max, k);
		swap_rows(B, index_max, k);
		for (size_t i = k+1; i < M; i++)
		{
			T mu = A[i][k] / A[k][k];
			for (size_t j = 0; j < N; j++)
			{
				A[i][j] = A[i][j] - mu*A[k][j];
			}
			B[i][0] = B[i][0] - mu*B[k][0];
		}
		
	}
	//*********direct order********
	//*********reverse order********
	X[M - 1][0] = B[M - 1][0] / A[M - 1][N - 1];
	for (int i = M-2; i >= 0; i--)
	{
		T summ = advanced_summ_of_row(A, X, i);
		X[i][0] = (B[i][0] - summ) / A[i][i];
	}
	//*********reverse order********
}

//********для своего класса матриц********//
template <typename T>
size_t find_max_of_column(Matrix<T> &A, const size_t &_j, const size_t &i_top)
{
	T max = fabs(A(i_top,_j)); 
	size_t i_max = i_top;
	for (size_t i = i_top; i < A.m; i++)
	{
		if (fabs(A(i, _j)) > max)
		{
			max = A(i, _j);
			i_max = i;
		}
	}
	return i_max;
}
template <typename T>
T advanced_summ_of_row(Matrix<T> &A, Matrix<T> &X, const size_t &_i)
{
	auto length = A.n;
	T summ = 0;
	for (size_t j = _i + 1; j < length; j++)
	{
		summ = summ + A(_i, j)*X(j, 0);
	}
	return summ;
}
template <typename T>
void method_Gauss_SLAU(Matrix<T> &A, Matrix<T> &B, Matrix<T> &X)
{
	//*********direct order********//
	const size_t m = A.m;
	const size_t n = A.n;
	for (size_t k = 0; k < m - 1; k++)
	{
		size_t i_max = find_max_of_column(A, k, k);
		A.swap_rows(i_max, k);
		B.swap_rows(i_max, k);
		for (size_t i = k + 1; i < m; i++)
		{
			float64 mu = A(i, k) / A(k, k);
			for (size_t j = 0; j < n; j++)
			{
				A(i,j) = A(i,j) - mu*A(k,j);
			}
			B(i,0) = B(i,0) - mu*B(k,0);
		}
	}
	//*********direct order********//
	//*********reverse order********//
	X(m - 1, 0) = B(m - 1, 0) / A(m - 1, n - 1);
	for (int i = m - 2; i >= 0; i--)
	{
		T summ = advanced_summ_of_row(A, X, i);
		X(i,0) = (B(i,0) - summ) / A(i,i);
	}
	//*********reverse order********//
}
//********для своего класса матриц********//