#pragma once
#include <Windows.h>
#include "Objective_function.h"
#ifndef _DEBUG
#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;
#endif

template <typename F,typename T>

class derivative_r //(f(x + h) - f(x)) / h
{
public:
	derivative_r(const F& f, const T& h) :f(f), h(h) {}
	Vertex2D<T> operator () (const T& x)
	{
		return (f(x + h) - f(x)) / h;
	}
	const Vertex2D<T> operator () (const T& x) const
	{
		return (f(x + h) - f(x)) / h;
	}
private:
	const F& f;
	T h;
};

template <typename F, typename T>

class derivative_l //(f(x) - f(x - h)) / h
{
public:
	derivative_l(const F& f, const T& h) :f(f), h(h) {}
	Vertex2D<T> operator () (const T& x)
	{
		return (f(x) - f(x - h)) / h;
	}
	const Vertex2D<T> operator () (const T& x) const
	{
		return (f(x) - f(x - h)) / h;
	}
private:
	const F& f;
	T h;
};

template <typename F, typename T>

class derivarive_2p //производная по двум точкам (слева и права) (f(x + h) - f(x - h)) / (2 * h)
{
public:
	derivarive_2p(const F& f, const T& h) :f(f), h(h) {}
	Vector2D<T> operator () (const T& x)
	{
		return (f(x + h) - f(x - h)) / (float(2) * h);
	}
	const Vector2D<T> operator () (const T& x) const
	{
		return (f(x + h) - f(x - h)) / (float(2) * h);
	}
private:
	const F& f;
	T h;
};

template<typename F, typename T>
class derivarive //производная по двум точкам (слева и права) (f(x + h) - f(x - h)) / (2 * h)
{
public:
	derivarive(const F& f, const T& h) :f(f), h(h) {}
	T operator () (vector<T>& var_arr, size_t i)
	{

		vector<T> delta_plus(var_arr);
		vector<T> delta_minus(var_arr);
		delta_plus[i] = delta_plus[i] + h;
		delta_minus[i] = delta_minus[i] - h;
		return (f(delta_plus) - f(delta_minus)) / (2*h);
/*		vector<T> delta_plus(var_arr);
		delta_plus[i] = delta_plus[i] + h;
		return ((f(delta_plus) - f(var_arr))) / h;*/
	}
	const T operator () (const vector<T>& var_arr, size_t i) const
	{
		/*vector<T> delta_plus(var_arr);
		vector<T> delta_minus(var_arr);
		delta_plus[i] = delta_plus[i] + h;
		delta_minus[i] = delta_minus[i] - h;
		return (f(delta_plus) - f(delta_minus)) / (2*h);*/
		vector<T> delta_plus(var_arr);
		delta_plus[i] = delta_plus[i] + h;
		return ((f(delta_plus) - f(var_arr))) / h;
	}
private:
	const F& f;
	T h;
};
template <typename F,typename T>
class second_derivative
{
public:
	second_derivative(F& f, const T& h) :h(h), fp(f, h), f(f) {}
	~second_derivative() {};
	T operator () (vector<T>& var_arr, size_t i, size_t di)
	{
		vector<T> delta_plus(var_arr);
		vector<T> delta_minus(var_arr);
		delta_plus[di] = delta_plus[di] + h;
		delta_minus[di] = delta_minus[di] - h;
		return (fp(delta_plus,i) - fp(delta_minus,i)) / (2 * h);
/*		vector<T> delta_plus(var_arr);
		delta_plus[di] = delta_plus[di] + h;
		return ((fp(delta_plus,i) - fp(var_arr,i))) / h;*/

	}
private:
	T h;
	F& f;
	derivarive<F, T> fp;
};

template <typename F,typename T>
//Упрощенное нахождение второй производной функции одной переменной
class second_derivative_simple
{
public:
	second_derivative_simple(F& f, const T& h) :h(h), f(f) {}
	~second_derivative_simple() {};
	T operator () (vector<T>& var_arr, size_t i, size_t di)
	{
		vector<T> delta_plus(var_arr);
		vector<T> delta_minus(var_arr);
		delta_plus[di] = delta_plus[di] + h;
		delta_minus[di] = delta_minus[di] - h;
		return (f(delta_plus) - 2 * f(var_arr) + f(delta_minus)) / (h*h);
	}
private:
	T h;
	F& f;
}; 

template <typename T>

class function
{
public:
	function() {};
	~function() {};

	T operator () (const vector<T> var_arr) const
	{
		T x1(var_arr[0]), x2(var_arr[1]);
		return (x1*x1 + 10*powf(x2 - sinf(x1),2));
	}
private:
};

template <typename T>

class temp_function
{
public:
	temp_function() {};
	~temp_function() {};

	T operator () (const vector<T> var_arr) const
	{
		T x1(var_arr[0]);
		float32 e = 2.71828;
		return (powf(x1,3) - x1 + powf(e,-x1));
	}
private:
};

void slau(double values[2][3], double& x1, double& x2)
{
	double delta = values[0][0] * values[1][1] - values[1][0] * values[0][1];
	double delta_x1 = values[1][1] * values[0][2] - values[0][1] * values[1][2];
	double delta_x2 = values[0][0] * values[1][2] - values[1][0] * values[0][2];
	x1 = delta_x1 / delta;
	x2 = delta_x2 / delta;
}
template<typename F, typename T>
void newton_minimization(F& f,vector<T>& variables )
{
	float32 h = 0.001;
	using d_f = derivarive<F, T>;
	using dd_f = second_derivative<F, T>;
	d_f dF(f, h);
	dd_f ddF(f, h);
	uint iterator = 0;
	float32& t = f.t;
	T& x1 = variables[0];
	T& x2 = variables[1];
	vector<T> var_new(2);
	T& x1n = var_new[0];
	T& x2n = var_new[1];
	t = f.curve.find_nearest(f.point);
	float32 alpha = 0.01;
	iterator++;
	auto g1 = dF(variables, 0);
	auto g2 = dF(variables, 1);
	auto G1 = ddF(variables, 0, 0);
	auto G2 = ddF(variables, 0, 1);
	auto G3 = ddF(variables, 1, 0);
	auto G4 = ddF(variables, 1, 1);
/*	auto g1 = 2 * x1 + 20 * (sinf(x1) - x2)*cosf(x1);
	auto g2 = 20 * (x2 - sinf(x1));
	auto G1 = 2 + 20 * (cosf(2 * x1) + x2*sinf(x1));
	auto G2 = -20 * cosf(x1);
	auto G3 = -20 * cosf(x1);
	auto G4 = 20;*/
	float64 p1, p2;

	float64 values[2][3] = {
		{ G1, G3, -g1 },
		{ G2, G4, -g2 }
	};
	slau(values, p1, p2);
	x1n = x1 + alpha * p1;
	x2n = x2 + alpha * p2;
	vector<T> residual; 
	vector<T> fv, fvn;
	fv.push_back(f(variables));
	fvn.push_back(f(var_new));
	vector<T> xAxsys;
	while (fabsf(f(var_new) - f(variables)) > 1e-6)
	{
		x1 = x1n;
		x2 = x2n;
		if ((iterator == 20)) alpha *= 10;
		t = f.curve.find_nearest(f.point);
		iterator++;
		g1 = dF(variables, 0);
		g2 = dF(variables, 1);
		G1 = ddF(variables, 0, 0);
		G2 = ddF(variables, 0, 1);
		G3 = ddF(variables, 1, 0);
		G4 = ddF(variables, 1, 1);
/*		g1 = 2 * x1 + 20 * (sinf(x1) - x2)*cosf(x1);
		g2 = 20 * (x2 - sinf(x1));
		G1 = 2 + 20 * (cosf(2 * x1) + x2*sinf(x1));
		G2 = -20 * cosf(x1);
		G3 = -20 * cosf(x1);
		G4 = 20;*/
		float64 values1[2][3] = {
			{ G1, G3, -g1 },
			{ G2, G4, -g2 }
		};
		slau(values1, p1, p2);
		if (x1 < 0) p1 *= -1;
		if (x2 < 0) p2 *= -1;
		x1n = x1 + alpha * p1;
		x2n = x2 + alpha * p2;
		residual.push_back(fabsf(f(var_new) - f(variables)));
		xAxsys.push_back(iterator);
		fv.push_back(f(variables));
		fvn.push_back(f(var_new));
#ifndef _DEBUG
		plt::clf();
		plt::plot(xAxsys, residual);
		plt::plot(xAxsys, fv);
		plt::plot(xAxsys, fvn);
		plt::pause(0.2);
#endif
		//Sleep(1000);
	}
}

template<typename T>
class function_bezier2point
{
public:
	function_bezier2point(BezierCurve<T>& _f,const  Vertex2D<T>& _point):f(_f),point(_point) {};
	~function_bezier2point() {};

	T operator () (const vector<T>& var) const
	{
		T t = var[0];
		return (f.getPoint(t).length(point));
	}

private:
	BezierCurve<T>& f;
	Vertex2D<T> point;
};
template <typename F, typename T>
void newton_minimization_simple(F& f, vector<T> &variables, const T left_border, const T right_border)
{
	float32 h(0.01);
	using d_f = derivarive<F, T>;
	using dd_f = second_derivative<F, T>;
	using dd_fs = second_derivative_simple<F, T>;
	d_f df(f, h);
	dd_f ddf(f, h);
	dd_fs ddfs(f, h);
	vector<T> new_variables;
	new_variables.push_back(0);
	T& x1 = variables[0];
	T& x1n = new_variables[0];
	auto z1 = df(variables, 0);
	auto z2 = ddf(variables, 0, 0);
	x1n = x1 - df(variables, 0) / ddf(variables, 0, 0);
	while (fabsf(f(new_variables) - f(variables)) > EPS)
	{
		x1 = x1n;
		x1n = x1 - df(variables, 0) / ddf(variables, 0, 0);
		if (x1 < left_border)	x1 = left_border;
		if (x1 > right_border)	x1 = right_border;
	}
	variables[0] = x1n;
}
template <typename F, typename T>
void method_bisection(F& f, vector<T> &variables, const T left_border, const T right_border)
{
	using d_f = derivarive<F, T>;
	T h = 0.001;
	d_f df(f, h);
	T& x = variables[0];
	T a = left_border;
	T b = right_border;
	x = (a + b) / 2;
	vector<float32> new_variables;
	new_variables.push_back(x);
	T& xn = new_variables[0];
	auto dfx = df(new_variables, 0);
	if (dfx < 0)
	{
		a = xn;
	}
	else
	{
		b = xn;
	}
	xn = (a + b) / 2;
	auto fx = f(variables);
	auto fxn = f(new_variables);
	while (fabsf(f(new_variables) - f(variables)) > (1e-10))
	{
		x = xn;
		dfx = df(new_variables, 0);
		if (dfx < 0)
		{
			a = xn;
		}
		else
		{
			b = xn;
		}
		xn = (a + b) / 2;
	}
	x = xn;
}
template <typename F, typename T>
void method_bisection_two_variables(F& f, vector<T> &variables, const T left_border, const T right_border)
{
	using d_f = derivarive<F, T>;
	T h = 0.001;
	d_f df(f, h);
	T& x = variables[0];
	variables[1] = x;
	T a = left_border;
	T b = right_border;
	x = (a + b) / 2;
	vector<float32> new_variables;
	new_variables.push_back(x);
	new_variables.push_back(x);
	T& xn = new_variables[0];
	auto dfx = df(new_variables, 0);
	if (dfx < 0)
	{
		a = xn;
	}
	else
	{
		b = xn;
	}
	xn = (a + b) / 2;
	new_variables[1] = xn;
	auto fx = f(variables);
	auto fxn = f(new_variables);
	while (fabsf(f(new_variables) - f(variables)) > (1e-6))
	{
		x = xn;
		variables[1] = x;
		dfx = df(new_variables, 0);
		if (dfx < 0)
		{
			a = xn;
		}
		else
		{
			b = xn;
		}
		xn = (a + b) / 2;
		new_variables[1] = xn;
	}
	x = xn;
}