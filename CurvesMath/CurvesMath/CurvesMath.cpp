// CurvesMath.cpp : Defines the entry point for the console application.
//
#include "stdafx.h"
#include "BezierCurve.h"
#include <iostream>
#include "MyMath.h"
#include <time.h>



int main()
{
	Vertex2D<float32> a(0,0);
	Vertex2D<float32> b(1, 1);
	Vertex2D<float32> c(2, 0);
	vector<Vertex2D<float32>> ppoints;
	ppoints.push_back(a);
	ppoints.push_back(b);
	ppoints.push_back(c);
	BezierCurve<float32> curve(ppoints);
	Vertex2D<float32> point1;
	Vertex2D<float32> dt1;
	//****************
	vector<Vertex2D<float32>> ppoints2;
	for (size_t i = 0; i < 9; i++)
	{
		Vertex2D<float32> X(float(3) - float(i)*4, float(i)*2 - float(3));
		ppoints2.push_back(X);
	}
	BezierCurve<float32> curve3(ppoints2);
	try
	{
		Vertex2D<float32> point = curve.getPoint(0);
		point = curve.getPoint(0.5);
		point1 = curve[0.5];
		dt1 = curve.dt(0.5);
		Vertex2D<float32> dot = curve.dt(0);
		dot = curve.dt(1);
		Vertex2D<float32> ddot = curve.ddt(0.5);
		ddot = curve.ddt(0);
		ddot = curve.ddt(1);
	}
	catch (const std::exception&)
	{
		cout << "Ошибка: неверный диапозон переменной" << endl;
	}
	BezierCurve<float32> curve2(ppoints);
	auto t = curve2.find_nearest(Vertex2D<float32>(1, 3));
	curve.increase();
	auto point = curve[0.5];
	auto dt = curve.dt(0.5);
	dt = curve.dt(1);
	//тут начинаем магию
	using d_bezier_f = derivative_r<BezierCurve<float32>, float32>;
	using dbf_l = derivative_l<BezierCurve<float32>, float32>;
	using dbf_2p = derivarive_2p<BezierCurve<float32>, float32>;
	using dd_bezier_f = derivative_r<d_bezier_f, float32>;
	using ddbf_l = derivative_l<dbf_l, float32>;
	using ddbf_2p = derivarive_2p<dbf_2p, float32>;
	auto h = 1e-3;
	//первый метод. справа
		d_bezier_f d_bezier_o(curve3, h);
		auto df = d_bezier_o(0.5);
		dd_bezier_f dd_bezier_o(d_bezier_o, h);
		auto ddf = dd_bezier_o(0.5);
	//второй метод слева
		dbf_l dbfo_l(curve3, h);
		auto df_l = dbfo_l(0.5);
		ddbf_l ddbfo_l(dbfo_l, h);
		auto ddf_l = ddbfo_l(0.5);
	//третий. двето точки по бокам
		dbf_2p dbfo_2p(curve3, h);
		auto df_2p = dbfo_2p(0.5);
		ddbf_2p ddbfo_2p(dbfo_2p, h);
		auto ddf_2p = ddbfo_2p(0.5);
	//аналитически
		auto DT = curve3.dt(0.5);
		auto DDT = curve3.ddt(0.5);
		auto XX = curve.curvature(0);
		cout << XX << endl;
		XX = curve.curvature(0.5);
		cout << XX << endl;
		XX = curve.curvature(0.6);
		cout << XX << endl;
		XX = curve.curvature(1);
		cout << XX << endl;
	//Vector2D test
		Vector2D<float32> vec3(3.0,3.0);
		vec3 = vec3 * 4;
		vec3 *= 0.5;
		vec3.normalize();
		auto l = vec3.length();
    return 0;
}

