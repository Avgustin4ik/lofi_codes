#pragma once
#include "Objective_function.h"

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

template <typename F, typename T>

class minimization_obj
{
public:
	minimization_obj() {};
	~minimization_obj() {};

	vector<T> Newton(BezierCurve<T>& _curve, const vector<T> argument_list,const Vertex2D<T> _point)
	{
		auto eps = 1e-5;
		auto alpha_1 = atan2f((_curve.dt(0.0).y) / (_curve.dt(0.0).x)) * 180 / PI;
		auto alpha_2 = fabs(atan2f((_curve.dt(1.0).y) / (_curve.dt(1.0).x)) * 180 / PI);//нужно следить, за получаемым числом
		auto gamma1 = argument_list[0];
		auto gamma2 = argument_list[1];
		auto m = _curve.m;
		const auto P0 = _curve.PPoints[0];
		const auto P3 = _curve.PPoints[3];
		//while (_curve.getPoint(near_point_arr[0]).length(_points[0] > eps) {}
		auto P1 = _curve.PPoints[0] + gamma1 * _curve.dt(0.0);
		auto P2 = _curve.PPoints[3] + gamma2 * _curve.dt(1.0);
		auto t = _curve.find_nearest(_point);
		Vertex2D<T> B = _curve.bernstein_data[0] * (powf(t, 0)*powf(1 - t, m - 0))*P0 + _curve.bernstein_data[1] * (powf(t, 0)*powf(1 - t, m - 1))*P1 + _curve.bernstein_data[2] * (powf(t, 0)*powf(1 - t, m - 2))*P2 + _curve.bernstein_data[3] * (powf(t, 0)*powf(1 - t, m - 3))*P3;
		B.length(_point);
	}
};
//для проверки


template<typename F, typename T>
class derivarive //производная по двум точкам (слева и права) (f(x + h) - f(x - h)) / (2 * h)
{
public:
	derivarive(const F& f, const T& h) :f(f), h(h) {}
	T operator () (vector<T>& var_arr, size_t i)
	{
/*		vector<T> delta_x(var_arr);
		delta_x[i] = delta_x[i] + h;
		return (f(delta_x) - f(var_arr)) / h;*/
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
	}
private:
	T h;
	F& f;
	derivarive<F, T> fp;
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
template<typename T>
void newton_minimization(objective_function<T>& f,vector<T>& variables )
{
	float32 h = 0.001;
	using d_f = derivarive<objective_function<float32>, float32>;
	using dd_f = second_derivative<objective_function<float32>, float32>;
	d_f dF(f, h);
	dd_f ddF(f, h);
	uint iterator = 0;
	float32& t = f.t;
	T& x1 = variables[0];
	T& x2 = variables[1];
	vector<T> variables_old;
	variables_old.push_back(1.0);
	variables_old.push_back(1.0);
	T& x1_old = variables_old[0];
	T& x2_old = variables_old[1];
	while (fabsf(f(variables) - f(variables_old)) > 1e-3)
	{
		x1_old = x1;
		x2_old = x2;
		t = f.curve.find_nearest(f.point);
		iterator++;
		auto g1 = dF(variables, 0);
		auto g2 = dF(variables, 1);
		auto G1 = ddF(variables, 0, 0);
		auto G2 = ddF(variables, 0, 1);
		auto G3 = ddF(variables, 1, 0);
		auto G4 = ddF(variables, 1, 1);
		float64 p1, p2;

		float64 values[2][3] = {
			{ G1, G3, -g1 },
			{ G2, G4, -g2 }
		};
		slau(values, p1, p2);
		x1 += p1;
		x2 += p2;
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
	float32 h(0.0001);
	using d_f = derivarive<F, T>;
	using dd_f = second_derivative<F, T>;
	d_f df(f, h);
	dd_f ddf(f, h);
	vector<T> new_variables;
	new_variables.push_back(0);
	T& x1 = variables[0];
	T& x1n = new_variables[0];
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