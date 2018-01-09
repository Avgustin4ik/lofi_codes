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
	
	
	vector < Vertex2D<float32>>	P;
	P.push_back(Vertex2D<float32>(0, 0));
	P.push_back(Vertex2D<float32>(1, 1));
	P.push_back(Vertex2D<float32>(2, 1));
	P.push_back(Vertex2D<float32>(3, 0));
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
	//***********************************
	vector<float32> variables;
	variables.push_back(0.5);
	variables.push_back(0.5);
	objective_function<float32> f(curve,point);
//	point = curve.getPoint(0.5);
//	newton_minimization<float32>(f, variables);

	vector<Vertex2D<float32>> P2(P);
	BezierCurve<float32> curve2(P);
	Vertex2D<float32> &p2 = curve2.PPoints[1];
	Vertex2D<float32> &p3 = curve2.PPoints[2];
	auto temp  = curve.dt(0);
	float32 gamma1 = sqrtf(2);
	p2 = curve2.PPoints[0] + temp * (gamma1*1.2);
	temp = curve.dt(1);
	temp.reverse();
	p3 = curve2.PPoints[3] + temp * (gamma1*1.2);
	vector<float32> Pnx, Pny;
	vector<float32> Bnx(101), Bny(101);
	for (size_t i = 0; i < 101; i++)
	{
		float32 t = float(i) / float(100);
		Bnx[i] = curve2.getPoint(t).x;
		Bny[i] = curve2.getPoint(t).y;
	}
	for (auto &i : curve2.PPoints)
	{
		Pnx.push_back(i.x);
		Pny.push_back(i.y);
	}

	vector<float32> near2X, near2Y;
	near2X.push_back(curve2.getPoint(curve2.find_nearest(point)).x);
	near2Y.push_back(curve2.getPoint(curve2.find_nearest(point)).y);
	//***********************************
	point = curve.getPoint(0.4);
	point = Vertex2D<float32>(1.25, 0.1);
	variables.clear();
	variables.push_back(1);
	variables.push_back(1);
	objective_function<float32> F(curve2, point);
	F.t = curve2.find_nearest(point);
	method_bisection_two_variables<objective_function<float32>, float32>(F, variables, 0, 2);
//	newton_minimization<objective_function<float32>, float32>(F, variables);
	gamma1 = variables[0];
	float32 gamma2 = variables[1];
	temp = curve.dt(0);
	p2 = curve2.PPoints[0] + temp * (gamma1);
	temp = curve.dt(1);
	temp.reverse();
	p3 = curve2.PPoints[3] + temp * (gamma2);
	Pnx.clear();	Pny.clear();	Bnx.clear();	Bny.clear();
	for (size_t i = 0; i < 101; i++)
	{
		float32 t = float(i) / float(100);
		Bnx.push_back(curve2.getPoint(t).x);
		Bny.push_back(curve2.getPoint(t).y);
	}
	for (auto &i : curve2.PPoints)
	{
		Pnx.push_back(i.x);
		Pny.push_back(i.y);
	}
/*	function<float32> primer;
	variables.clear();
	variables.push_back(1);
	variables.push_back(1);
	newton_minimization<function<float32>, float32>(primer, variables);*/

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
	plt::plot(near2X, near2Y, "X-y");
	plt::grid(true);
	plt::axis("equal");
	plt::show();
#endif
	
	

    return 0;
	
}

