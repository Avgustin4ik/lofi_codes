#pragma once
#include "Objective_function.h"
#include "Gauss_math.h"
#ifndef _DEBUG
#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;
#endif
template<typename F, typename T>
class derivarive //производная по двум точкам (слева и права) (f(x + h) - f(x - h)) / (2 * h)
{
public:
	derivarive(F& f, const T& h) :f(f), h(h) {}
	T operator () (vector<T>& var_arr, size_t i)
	{

		vector<T> delta_plus(var_arr);
		vector<T> delta_minus(var_arr);
		delta_plus[i] = delta_plus[i] + h;
		delta_minus[i] = delta_minus[i] - h;
		return (f(delta_plus) - f(delta_minus)) / (2*h);
	}
	const T operator () (const vector<T>& var_arr, size_t i) const
	{
		vector<T> delta_plus(var_arr);
		vector<T> delta_minus(var_arr);
		delta_plus[i] = delta_plus[i] + h;
		delta_minus[i] = delta_minus[i] - h;
		return (f(delta_plus) - f(delta_minus)) / (2*h);
	}

	Matrix<T> operator ()(const vector<T>& var)
	{
		size_t size = var.size();
		Matrix<T> result(size, 1);
		for (size_t i = 0; i < size; i++)
		{
			vector<T> delta_plus(var);
			vector<T> delta_minus(var);
			delta_plus[i] = delta_plus[i] + h;
			delta_minus[i] = delta_minus[i] - h;
			T a = (f(delta_plus) - f(delta_minus)) / (2 * h);
			result(i, 0) = EQUAL(a, 0.0);
		}
		return result;
	}

private:
	F& f;
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
	}
	Matrix<T> operator ()(const vector<T>& var)
	{
		size_t size = var.size();
		Matrix<T> result(size, size);
		for (size_t i = 0; i < size; i++)
			for (size_t j = 0; j < size; j++)
			{
				size_t index = i*size + j;
				vector<T> delta_plus(var);
				vector<T> delta_minus(var);
				delta_plus[j] = delta_plus[j] + h;
				delta_minus[j] = delta_minus[j] - h;
				T a = (fp(delta_plus)(i,0) - fp(delta_minus)(i,0)) / (2 * h);
				result(i, j) = EQUAL(a,0.0);
			}
		return result;
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
		return (x1*x1*5+x1*3);
	}
private:
};

template <typename T>

void slau(T values[2][3], T& x1, T& x2)
{
	T delta = values[0][0] * values[1][1] - values[1][0] * values[0][1];
	T delta_x1 = values[1][1] * values[0][2] - values[0][1] * values[1][2];
	T delta_x2 = values[0][0] * values[1][2] - values[1][0] * values[0][2];
	x1 = delta_x1 / delta;
	x2 = delta_x2 / delta;
}

template<typename F, typename T>

void newton_minimization(objective_function_tangent<T>& f,vector<T>& variables, Configuration _config )
{
	vector<T> initial_data(variables);
	vector<T> npx, npy, fx, fy, dfx, dfy;
	vector<float32>  x;
	vector<float32>  y;
	x.push_back(f.point.x);
	y.push_back(f.point.y);
	T h = _config.h;
	float32 alpha = _config.alpha;
	bool isSolutionNotReached = false;
	derivarive<F, T> df(f,h);
	second_derivative<F, T> ddf(f, h);
	vector<T> variables_new(2);
	Matrix<T> p(variables.size(), 1);
	
	vector<float32> Pnx, Pny, Bnx, Bny;//для отображения
#ifndef _DEBUG
	plt::clf();
	for (auto i = 0; i < 101; i++) {
		auto z = float(i) / float(100);
		Bnx.push_back(f.curve.getPoint(z).x);
		Bny.push_back(f.curve.getPoint(z).y);
	}
	for (auto &i : f.curve.PPoints) {
		Pnx.push_back(i.x);
		Pny.push_back(i.y);
	}
	plt::grid(true);
	plt::axis("equal");
	plt::plot(x, y, "D");
	plt::plot(Pnx, Pny, "x--r");
	plt::plot(Bnx, Bny);
	plt::pause(0.5);
#endif
	f.recompute(variables);
	auto f_v = f(variables);
	size_t size = variables.size();
	auto g(df(variables));
	auto G(ddf(variables));
	g = g * -1;
	method_Gauss_SLAU(G, g, p);
	for (size_t i = 0; i < p.m; i++)
	{
		size_t j = 0;
		variables_new[i] = variables[i] + alpha * p(i, j);
	}
	auto f_vn = f(variables_new);
	f.recompute(variables_new);
	int i = 1;
	int step = 0;
	Pnx.clear(); Pny.clear(); Bnx.clear(); Bny.clear();
	bool norm(true);
#ifndef _DEBUG
	plt::clf();
	for (auto i = 0; i < 101; i++) {
		auto z = float(i) / float(100);
		Bnx.push_back(f.curve.getPoint(z).x);
		Bny.push_back(f.curve.getPoint(z).y);
	}
	for (auto &i : f.curve.PPoints) {
		Pnx.push_back(i.x);
		Pny.push_back(i.y);
	}
	plt::grid(true);
	plt::axis("equal");
	plt::plot(x, y, "D");
	plt::plot(Pnx, Pny, "x--r");
	plt::plot(Bnx, Bny);
	plt::pause(0.5);
#endif
	vector<T> xx, yy;
	while (fabs(f(variables_new)) > 1e-5)
	//while (powf((f.getBezierPoint(variables_new).x - f.point.x), 2) + powf((f.getBezierPoint(variables_new).y - f.point.y), 2) > 1e-5)
	{
		T &x1n = variables_new[0];
		T &x2n = variables_new[1];
		size = variables.size();
		Matrix<T> p(size, 1);
		vector<Vertex2D<T>> &PP = f.curve.PPoints;
		i++;
		int m = f.curve.PPoints.size()-1;
		if (i > (m-2)*40) { 
			variables = initial_data;
			f.recompute(variables);
			f.add_PPoint(variables);
			initial_data = variables;
			p = Matrix<T>(variables.size(), 1);
			step = 0;
			alpha = _config.alpha;
			variables_new.resize(size + 2);
			norm = false;
		}
		if (i == _config.iterations_limit) { isSolutionNotReached = true; break; }
		if (PP[1].x >= PP[2].x || PP[PP.size() - 2].x <= PP[PP.size() - 3].x) {
			f.recompute(variables);
			f.add_PPoint(variables);
			initial_data = variables;
			p = Matrix<T>(variables.size(), 1);
			step = 0;
			alpha = _config.alpha;
			variables_new.resize(size + 2);
			norm = false;
		}
		if (((x1n) <= 0.1) || ((x2n) <= 0.1)) {
			f.recompute(variables);
			f.add_PPoint(variables);
			initial_data = variables;
			p = Matrix<T>(variables.size(), 1);
			step = 0;
			alpha = _config.alpha;
			variables_new.resize(size + 2);
			norm = false;
		} 
		step++;
		if (norm) variables = variables_new;
		norm = true;
		if (step > (m - 2) * 30) alpha = 1;
		auto g(df(variables));
		auto G(ddf(variables));
		g = g * -1;
		method_Gauss_SLAU(G, g, p);
		for (size_t i = 0; i < p.m; i++)
		{
			variables_new[i] = variables[i] + alpha * p(i, 0);
		}
		f_vn = f(variables_new);
		f.recompute(variables_new);
		fx.push_back(i);
		fy.push_back(f_vn);
		Pnx.clear(); Pny.clear(); Bnx.clear(); Bny.clear(); npx.clear(); npy.clear(); dfx.clear(); dfy.clear();
#ifndef _DEBUG
		plt::clf();
		plt::subplot(2, 3, 1);
		for (auto i = 0; i < 101; i++) {
			auto z = float(i) / float(100);
			auto P = f.curve.getPoint(z);
			Bnx.push_back(P.x);
			Bny.push_back(P.y);
		}
		for (auto &i : f.curve.PPoints) {
			Pnx.push_back(i.x);
			Pny.push_back(i.y);
		}
		npx.push_back(f.curve.getPoint(f.curve.find_nearest(f.point)).x);
		npy.push_back(f.curve.getPoint(f.curve.find_nearest(f.point)).y);
		plt::plot(npx, npy,"D");
		plt::grid(true);
		plt::axis("equal");
		plt::plot(x, y, "D");
		plt::plot(Pnx, Pny, "x--r");
		plt::plot(Bnx, Bny);

		plt::subplot(2, 3, 2);
		plt::plot(fx, fy);
		plt::grid(true);
		
		plt::subplot(2, 3, 3);
		plt::grid(true);
		plt::axis("equal");
		auto c = f.curve;
		dfx.push_back(0);
		dfx.push_back(c.dt(0).x);
		dfy.push_back(0);
		dfy.push_back(c.dt(c.find_nearest(f.point)).y);
		plt::plot(dfx, dfy);
		plt::grid(true);
		plt::axis("equal");
		plt::pause(0.001);
#endif
	}
	if (!isSolutionNotReached)
	{
		variables = variables_new;
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
	while (fabsf(f(new_variables) - f(variables)) > (1e-8))
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
	f.recompute(new_variables[0],new_variables[1]);
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
		f.recompute(new_variables[0], new_variables[1]);
		//****************************для отображения
		vector<float32> Pnx, Pny, Bnx, Bny;
		for (size_t i = 0; i < 101; i++)
		{
			float32 t = float(i) / float(100);
			Bnx.push_back(f.curve.getPoint(t).x);
			Bny.push_back(f.curve.getPoint(t).y);
		}
		for (auto &i : f.curve.PPoints)
		{
			Pnx.push_back(i.x);
			Pny.push_back(i.y);
		}
#ifndef _DEBUG
		plt::clf();
		//		plt::plot(xAxsys, residual);
		//		plt::plot(xAxsys, fv);
		//		plt::plot(xAxsys, fvn);
		plt::grid(true);
		plt::axis("equal");
		vector<float32> x, y;
		x.push_back(f.point.x);
		y.push_back(f.point.y);
		plt::plot(x, y, "D");
		plt::plot(Pnx, Pny, "x--r");
		plt::plot(Bnx, Bny);
		plt::pause(0.5);
#endif
		//****************************для отображения
		
	}
	x = xn;
	variables[1] = x;
}


