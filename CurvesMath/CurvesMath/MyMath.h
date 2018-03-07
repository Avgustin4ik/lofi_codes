#pragma once
#include "Objective_function.h"
#include "Gauss_math.h"
#include "Circle.h"
#include "FishBone.h"
#include "SidesFunction.h"
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
		return (f(delta_plus) - f(delta_minus)) / (2 * h);
	}
	const T operator () (const vector<T>& var_arr, size_t i) const
	{
		vector<T> delta_plus(var_arr);
		vector<T> delta_minus(var_arr);
		delta_plus[i] = delta_plus[i] + h;
		delta_minus[i] = delta_minus[i] - h;
		return (f(delta_plus) - f(delta_minus)) / (2 * h);
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


template <typename F, typename T>
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
		return (fp(delta_plus, i) - fp(delta_minus, i)) / (2 * h);
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
				T a = (fp(delta_plus)(i, 0) - fp(delta_minus)(i, 0)) / (2 * h);
				result(i, j) = EQUAL(a, 0.0);
			}
		return result;
	}
private:
	T h;
	F& f;
	derivarive<F, T> fp;
};

template <typename F, typename T>
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
void minimization_Newthon(of_CamberLine<T>& f, vector<T>& variables, Configuration _config)
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
	derivarive<F, T> df(f, h);
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
	obj_function_alpha<of_CamberLine<T>, T> f_alpha(f, variables, p_vec);
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
	plt::pause(0.01);
#endif
	vector<T> xx, yy;
	while (fabs(f(variables_new)) > 1e-3)
		//while (powf((f.getBezierPoint(variables_new).x - f.point.x), 2) + powf((f.getBezierPoint(variables_new).y - f.point.y), 2) > 1e-5)
	{
		size = variables.size();
		Matrix<T> p(size, 1);
		vector<Vertex2D<T>> &PP = f.curve.PPoints;
		i++;
		int m = f.curve.PPoints.size() - 1;
		if (i > (m - 2) * 30) {
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
		if ((powf(variables_new[0], 2) <= 0.001) || (powf(variables_new[1], 2) <= 0.001)) {
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
		plt::plot(npx, npy, "D");
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

template<typename F, typename T>
void minimization_Newthon(F &f, vector<T> &x, Configuration _config)
{
	vector<double> fx, fy;
	int i = 0;
	auto initialData(x);
	vector<T> alpha;
	alpha.push_back(0.1);
	while (fabsf(f(x)) > 1e-3)
	{
		vector<Vertex2D<float32>> &cP(f.getCurve().PPoints);
		if (i % 50 == 0 && i != 0) (f.increaseCurve(x));
		if (i == _config.iterations_limit) break;
		if (cP[0].length(cP[1])<1e-3 || cP[cP.size() - 1].length(cP[cP.size() - 2])<1e-3)
		{
			x = initialData;
			f.increaseCurve(x);
			initialData = x;
		}
		for (auto &ptr : x)
			if (ptr > 1)
			{
				x = initialData;
				f.increaseCurve(x);
				initialData = x;
			}
		T h = _config.h;
		derivarive<F, T> df(f, h);
		second_derivative<F, T> ddf(f, h);
		Matrix<T> p(x.size(), 1);
		Matrix<T> g(df(x));
		Matrix<T> G(ddf(x));
		matrixScaling(G, g);
		method_Gauss_SLAU(G, g, p);

		vector<T> p_vec;
		p_vec.reserve(p.m);
		for (size_t i = 0; i < p.m; i++)
			p_vec.emplace_back(p(i, 0));
		obj_function_alpha<SidesFunction<float32>, float32> f_alpha(f, x, p_vec);
		method_bisection<obj_function_alpha<SidesFunction<float32>, float32>, T>(f_alpha, alpha, 0.0, 1.0);
		for (size_t i = 0; i < p.m; i++)
		{
			size_t j = 0;
			x[i] = x[i] + alpha[0] * p(i, j);
		}
		i++;
		vector<double> Pnx, Pny, Bnx, Bny, cpx, cpy;
#ifndef _DEBUG
		plt::clf();
		plt::subplot(2, 1, 1);
		
		plot_curve(f.getCurve());
		plot_fishbones(f.getFishBones());
		for (auto &i : f.getConditionPoints())
		{
			cpx.push_back(i.x);
			cpy.push_back(i.y);
		}
		plt::grid(true);
		plt::axis("equal");
		plt::plot(Pnx, Pny, "x--r");
		plt::plot(Bnx, Bny);
		plt::plot(cpx, cpy, "D");
		plt::subplot(2, 1, 2);
		plt::grid(true);
		fx.push_back(i);
		fy.push_back(f(x));
		plt::plot(fx, fy, "-b");
		plt::pause(0.001);
#endif // !_DEBUG

	}

}

#ifndef _DEBUG
template<typename F, typename T>
void minimize_cpp(BfgsSolver<SidesFunction<float32>> &solver, F &f, VectorXd &x)
{
	VectorXd initialData(x);
	vector<T> variables(-1, x.size());
	while (f.value(x) > 1e-3)
	{
		solver.minimize(f, x);
		plot_curve(f.getCurve());
		plt::pause(0.1);
		std::cout << "argmin      " << x.transpose() << std::endl;
		if (f.value(x) > 1e-3)
		{
			x = initialData;
			for (size_t i = 0; i < x.size(); i++)	variables[i] = x[i];
			f.increaseCurve(variables);
			cout << "Point is added" << endl;
			x.resize(variables.size());
			for (size_t i = 0; i < x.size(); i++)	x[i] = variables[i];
			std::cout << "argmin after adding point     " << x.transpose() << std::endl;
			initialData = x;
		}
	}
}

#endif // !_DEBUG

template<typename T>
class function_bezier2point
{
public:
	function_bezier2point(BezierCurve<T>& _f, const  Vertex2D<T>& _point) :f(_f), point(_point) {};
	~function_bezier2point() {};

	T operator () (const vector<T>& var) const
	{
		T t = var[0];
		t = EQUAL2EPS(t, 0.0, 1e-3);
		t = EQUAL2EPS(t, 1.0, 1e-3);
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
	f.recompute(new_variables[0], new_variables[1]);
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

	while (f(variables) > 1e-3)
	{
		index = 0;
		iterations++;
		while (index < variables.size())
		{
			method_bisection<obj_functionCoorDescent<float32>, float32>(f_obj, variables, 0.0, 1.0);
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
		if (iterations > (P.size() - 3) * 50)
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
	if (omega > 180 || omega < 0) throw(exception());
	Vector2D<T> tangent;
	Vertex2D<T> P;
	vector<Vertex2D<T>> points;
	points.reserve(2);
	float64 _omega(omega);
	if (!isLeadingEdge) _omega *= -1;
	float64 angle = PI / 2 - RAD(_omega / 2);
	if (isLeadingEdge) {
		P = curve.PPoints[0];
		tangent = curve.dt(0);
	}
	else {
		P = curve.PPoints[curve.PPoints.size() - 1];
		tangent = curve.dt(1);
	}
	float64 phi = PI - angle;
	Vector2D<T> direct1 = tangent + phi;
	Vector2D<T> direct2 = tangent - phi;
	points.emplace_back(P + direct1 * radius);
	points.emplace_back(P + direct2 * radius);
	return points;
}


//#ifndef _DEBUG
//using namespace cppoptlib;
//using Eigen::VectorXd;
//template<typename F, typename S>
//void minimization(F &f, VectorXd &x,ISolver<f,0> &solver)
//{
//	while (f(x)>1e-3)
//	{
//		solver.minimize(f, x);
//
//
//	}
//}
//
//
//
//
//
//
//#endif // !_DEBUG

