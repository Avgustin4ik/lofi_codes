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
		for (size_t j = 0; j < size; j++)
			for (size_t i = 0; i < size; i++)
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

template<typename F, typename T>
void minimization_Newthon(of_CamberLine<T>& f,vector<T>& variables, Configuration _config )
{
	vector<T> initial_data(variables);
	vector<T> npx, npy, fx, fy, dfx, dfy, ax, ay, k1y, k2y, curvX, curvY;
	vector<float32>  x;
	vector<float32>  y;
	x.push_back(f.point.x);
	y.push_back(f.point.y);
	T h = _config.h;
	vector<float32> alpha(1, _config.alpha);
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
	Matrix<T> g(df(variables));
	Matrix<T> G(ddf(variables));
	matrixScaling(G, g);
	g = g * -1;
	method_Gauss_SLAU(G, g, p);

	vector<T> p_vec;
	p_vec.reserve(p.m);
	for (size_t i = 0; i < p.m; i++)
		p_vec.emplace_back(p(i, 0));
	obj_function_alpha<T> f_alpha(f, variables, p_vec);
	method_bisection(f_alpha, alpha, LEFT_BORDER, RIGHT_BORDER);

	for (size_t i = 0; i < p.m; i++)
	{
		size_t j = 0;
		variables_new[i] = variables[i] + alpha[0] * p(i, j);
	}
	auto f_vn = f(variables_new);
	f.recompute(variables_new);
	int i = 1;
	Pnx.clear(); Pny.clear(); Bnx.clear(); Bny.clear(); ax.clear(); ay.clear();
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
	while (fabs(f(variables_new)) > 1e-3)
	//while (powf((f.getBezierPoint(variables_new).x - f.point.x), 2) + powf((f.getBezierPoint(variables_new).y - f.point.y), 2) > 1e-5)
	{
		size = variables.size();
		Matrix<T> p(size, 1);
		vector<Vertex2D<T>> &PP = f.curve.PPoints;
		i++;
		int m = f.curve.PPoints.size()-1;
		if (i > (m-2)*50) { 
			variables = initial_data;
			f.recompute(variables);
			f.add_PPoint(variables);
			initial_data = variables;
			p = Matrix<T>(variables.size(), 1);
			variables_new.resize(variables.size());
			norm = false;
		}
		if (i == _config.iterations_limit) { isSolutionNotReached = true; break; }
		if (PP[1].x >= PP[2].x || PP[PP.size() - 2].x <= PP[PP.size() - 3].x) {
			f.recompute(variables);
			f.add_PPoint(variables);
			initial_data = variables;
			p = Matrix<T>(variables.size(), 1);
			variables_new.resize(variables.size());
			norm = false;
		}
		if ((powf(variables_new[0],2) <= 0.001) || (powf(variables_new[1],2) <= 0.001)) {
			variables = initial_data;
			f.recompute(variables);
			f.add_PPoint(variables);
			initial_data = variables;
			p = Matrix<T>(variables.size(), 1);
			variables_new.resize(variables.size());
			norm = false;
		} 
		if (norm) variables = variables_new;
		norm = true;
		Matrix<T> g(df(variables));
		Matrix<T> G(ddf(variables));
		auto _G = G.invers();
		auto DET = G.determinant();
		matrixScaling(G, g);
		g = g * -1;
		method_Gauss_SLAU(G, g, p);

		p_vec.clear();
		p_vec.reserve(p.m);
		for (size_t i = 0; i < p.m; i++)
			p_vec.emplace_back(p(i, 0));
		obj_function_alpha<T> f_alpha(f, variables, p_vec);
		method_bisection(f_alpha, alpha, LEFT_BORDER, RIGHT_BORDER);

		for (size_t i = 0; i < p.m; i++)
		{
			variables_new[i] = variables[i] + alpha[0] * p(i, 0);
		}
		f_vn = f(variables_new);
		f.recompute(variables_new);
		fx.push_back(i);
		fy.push_back(f_vn);
		Pnx.clear(); Pny.clear(); Bnx.clear(); Bny.clear(); npx.clear(); npy.clear(); dfx.clear(); dfy.clear(), curvX.clear(), curvY.clear();
#ifndef _DEBUG
		plt::clf();
		plt::subplot(3, 2, 1);
		for (auto i = 0; i < 101; i++) {
			auto z = float(i) / float(100);
			auto P = f.curve.getPoint(z);
			Bnx.push_back(P.x);
			Bny.push_back(P.y);
			curvX.push_back(z);
			curvY.push_back(f.curve.curvature(z));
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

		plt::subplot(3, 2, 2);
		plt::plot(fx, fy);
		plt::grid(true);
		
		plt::subplot(3, 2, 3);
		plt::grid(true);
		plt::axis("equal");
		auto c = f.curve;
		dfx.push_back(0);
		dfx.push_back(c.dt(0).x);
		dfy.push_back(0);
		dfy.push_back(c.dt(c.find_nearest(f.point)).y);
		plt::plot(dfx, dfy);
		plt::grid(true);

		plt::subplot(3, 2, 4);
		plt::ylim(0, 1);
		plt::ylabel("alpha");
		ax.push_back(i);
		ay.push_back(alpha[0]);
		plt::plot(ax, ay);
		plt::grid(true);

		plt::subplot(3, 2, 5);
		plt::ylim(-5, 1);
		plt::ylabel("curvature");
		k1y.push_back(f.k1);
		k2y.push_back(f.k2);
		plt::plot(ax, k1y);
		plt::plot(ax, k2y);
		plt::grid(true);

		plt::subplot(3, 2, 6);
		plt::ylabel("curvature");
		plt::plot(vector<double>({ 0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0 }), f.curvature);
		plt::plot(curvX, curvY);
		
		plt::grid(true);

		plt::pause(0.1);
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
void minimization_NewthonSimple(F& f, vector<T> &variables, const T left_border, const T right_border)
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
	while (fabsf(f(new_variables) - f(variables)) > (1e-4))
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
template <typename F, typename T>

void method_bisection(obj_functionCoorDescent<T> &f, vector<T> &variables, const T left_border, const T right_border)
{
	using d_f = derivarive<obj_functionCoorDescent<T>, T>;
	T h = 0.001;
	d_f df(f, h);
	int index = f.index;
	T& x = variables[index];
	T a = left_border;
	T b = right_border;
	x = (a + b) / 2;
	vector<float32> new_variables(variables);
	new_variables[index] = x;
	T& xn = new_variables[index];
	auto dfx = df(new_variables, index);
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
	while (fabsf(f(new_variables) - f(variables)) > (1e-5))
	{
		x = xn;
		dfx = df(new_variables, index);
		if (dfx < 0)
		{
			a = xn;
		}
		else
		{
			b = xn;
		}
		xn = (a + b) / 2;
		f.f.recompute(new_variables);
	}
	x = xn;
}


template <typename F, typename T>
void minimization_CoordinateDescent(of_CamberLine<T>& f, vector<T> &variables, Configuration _config)
{
	vector<float32> Bnx, Bny, Pnx, Pny, fx, fy, npx, npy;
	vector<float32> x(1, f.point.x), y(1, f.point.y);
	derivarive<F, T> df(f, _config.h);
	int index = 0;
	obj_functionCoorDescent<T> f_obj(f, index);
	vector<T> initial_data(variables);
	int iterations(0);
	bool isSolution = true;
	vector<Vertex2D<T>> &P = f.curve.PPoints;
	
	while (f(variables)>1e-3)
	{
		index = 0;
		iterations++;
		while (index < variables.size())
		{
			method_bisection<obj_functionCoorDescent<float32>, float32>(f_obj, variables, 0.0,1.0);
			index++;
			f.recompute(variables);
		}
		if (iterations == _config.iterations_limit) { isSolution = false; break; }
		if ((P[0].x >= P[1].y) || (P[P.size() - 3].x < P[P.size() - 4].x))
		{
			variables = initial_data;
			f.recompute(variables);
			f.add_PPoint(variables);
			initial_data = variables;
		}
		if ((variables[0] < 0.005) || (variables[variables.size() - 1] < 0.005))
		{
			variables = initial_data;
			f.recompute(variables);
			f.add_PPoint(variables);
			initial_data = variables;
		}
		if (iterations > (P.size() - 3)*50) 
		{
			variables = initial_data;
			f.recompute(variables);
			f.add_PPoint(variables);
			initial_data = variables;
		}
		Bnx.clear(); Bny.clear(); Pnx.clear(); Pny.clear(); npx.clear(); npy.clear();
#ifndef _DEBUG
		plt::clf();
		plt::subplot(2, 2, 1);
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
		plt::plot(npx, npy, "D");
		plt::grid(true);
		plt::axis("equal");
		plt::plot(x, y, "D");
		plt::plot(Pnx, Pny, "x--r");
		plt::plot(Bnx, Bny);
		
		plt::subplot(2, 2, 2);
		fx.push_back(iterations); fy.push_back(f(variables));
		plt::plot(fx, fy);
		plt::grid(true);

		plt::pause(0.01);
#endif
	}
	
	
}
template <typename T>
vector<Vertex2D<T>> getEdgePoints(const bool &isLeadingEdge, BezierCurve<T> &curve, const T &omega, const T &radius)
{
	
	Vector2D<T> tangent;
	Vertex2D<T> P;
	vector<Vertex2D<T>> points;
	points.reserve(2);
	float64 angle = 90 - omega / 2;
	angle = angle * 180 / PI;
	if (isLeadingEdge) {
		P = curve.PPoints[0];
		tangent = curve.dt(0);
	}
	else {
		P = curve.PPoints[curve.PPoints.size() - 1];
		tangent = curve.dt(1);
	}
	float64 a = atan(tangent.y / tangent.x);
	float64 phi = 2 * PI + a - angle;
	Vector2D<T> direct1 = tangent + phi;
	Vector2D<T> direct2 = tangent - phi;
	points.emplace_back(P + direct1 * radius);
	points.emplace_back(P + direct2 * radius);
	return points;
}

