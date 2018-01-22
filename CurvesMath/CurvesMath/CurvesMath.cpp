// CurvesMath.cpp : Defines the entry point for the console application.
//
#include "stdafx.h"
#include "BezierCurve.h"
#include <iostream>
#include "MyMath.h"
#include <time.h>

#ifndef _DEBUG
#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;
#endif

int main()
{
	Configuration config;
	vector < Vertex2D<float32>>	P;
	P.push_back(Vertex2D<float32>(0, 0));
	P.push_back(Vertex2D<float32>(0.5 * cos(60*PI/180.0), 0.5 *sin(60*PI/180.0)));
	P.push_back(Vertex2D<float32>(1 - 0.5 *cos(21.5 * PI / 180.0), 0.5 *sin(60 * PI / 180.0)));
	P.push_back(Vertex2D<float32>(1, 0));
	BezierCurve<float32> curve(P);
	std::vector<float> x_(10);
	std::vector<float> y_(10);
	vector<float32> Bx(101), By(101);
	Vertex2D<float32> point(1, 0.5);
	for (size_t i = 0; i < 101; i++)
	{
		float32 t = float(i) / float(100);
		Bx[i] = curve.getPoint(t).x;
		By[i] = curve.getPoint(t).y;
	}
	vector<float32> Px;
	vector<float32> Py;
	for (auto &i : curve.PPoints)
	{
		Px.push_back(i.x);
		Py.push_back(i.y);
	}
	Vertex2D<float32> nearPoint = curve.getPoint(curve.find_nearest(point));
	vector<float32> nearX, nearY;
	nearX.push_back(nearPoint.x);
	nearY.push_back(nearPoint.y);

	

	using POINTS = vector<Vertex2D<float32>>;
	using PP = Vertex2D<float32>;
	BezierCurve<float32> line_curve(POINTS({PP(0.0,0.0), PP(1.0,0.0)}));
	point = Vertex2D<float32>(0.16, 0.11);
	vector<float32> variables;
	variables.push_back(0.5);
	variables.push_back(0.5);
	objective_function_tangent<float32> F(curve,line_curve, point);


	newton_minimization<objective_function_tangent<float32>, float32>(F, variables, config);

	/*float32 g1 = variables[0];
	float32 g2 = variables[1];
	auto &p2 = curve.PPoints[1];
	auto &p3 = curve.PPoints[2];
	auto vec = curve.dt(0);
	p2 = curve.PPoints[0] + vec*powf(g1,2);
	vec = curve.dt(1);
	vec.reverse();
	p3 = curve.PPoints[3] + vec*powf(g2,2);*/
	vector<float32> Pnx, Pny, Bnx, Bny;
	for (size_t i = 0; i < 101; i++)
	{
		float32 t = float(i) / float(100);
		Bnx.push_back(curve.getPoint(t).x);
		Bny.push_back(curve.getPoint(t).y);
	}
	for (auto &i : curve.PPoints)
	{
		Pnx.push_back(i.x);
		Pny.push_back(i.y);
	}
#ifndef _DEBUG
	plt::figure();
	plt::plot(Px, Py, "D--g");
	plt::plot(Bx, By);
	plt::plot(Pnx, Pny, "x--r");
	plt::plot(Bnx, Bny);
	vector<float32> x, y;
	x.push_back(point.x);
	y.push_back(point.y);
	plt::plot(x,y,"D");
	plt::plot(nearX, nearY, "X-r");
//	plt::plot(near2X, near2Y, "X-y");
	plt::grid(true);
	plt::axis("equal");
	plt::show();

#endif
	

    return 0;
	
}

