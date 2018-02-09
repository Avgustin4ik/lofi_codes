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

using namespace cppoptlib;
using Eigen::VectorXd;


int main()
{
	
	
	Configuration config;
	vector < Vertex2D<float32>>	P;
	P.push_back(Vertex2D<float32>(0, 0));
	P.push_back(Vertex2D<float32>(0.5 * cos(60*PI/180.0), 0.5 *sin(60*PI/180.0)));
	P.push_back(Vertex2D<float32>(1 - 0.5 *cos(21.5 * PI / 180.0), 0.5 *sin(21.5 * PI / 180.0)));
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
	
	//point = Vertex2D<float32>(0.4, 0.11);
	vector<float32> variables;
	variables.push_back(0.5);
	variables.push_back(0.5);
	of_CamberLine<float32> f(curve,line_curve,point);
	CamberLineFunction<float32> fun(curve, line_curve, point);
	//minimization_Newthon<of_CamberLine<float32>, float32>(f, variables, config);
	//minimization_CoordinateDescent<of_CamberLine<float32>, float32>(F, variables, config);

	
	BfgsSolver<CamberLineFunction<float32>> solver;
	//LbfgsbSolver<CamberLineFunction<float32>> solver;
	//LbfgsSolver<CamberLineFunction<float32>> solver;
	VectorXd X(2); X << 0.5, 0.5;
	Criteria<double> stopCriteria;
	stopCriteria.iterations = 50;
	stopCriteria.xDelta = 1e-4;
	stopCriteria.fDelta = 1e-4;
	solver.setStopCriteria(stopCriteria);

	VectorXd initialData(X);
	while (fun(X) > 1e-3)
	{
		std::cout << "argmin befor     " << X.transpose() << std::endl;
		solver.minimize(fun, X);
		for (size_t i = 0; i < X.size(); i++)
		{
			variables[i] = X[i];
		}
		fun.recompute(variables);
		cout << solver.status() << endl;
		solver.criteria().print(cout);
		std::cout << "argmin after     " << X.transpose() << std::endl;
		std::cout << "f in argmin " << fun(X) << std::endl;
		if (fun(X) > 1e-3) {
			X = initialData;
			fun.add_PPoint(X);
			cout << "Point is added" << endl;
			std::cout << "argmin after adding point     " << X.transpose() << std::endl;
			initialData = X;
			variables.resize(X.size());
		}
	}
	for (size_t i = 0; i < X.size(); i++)
	{
		variables[i] = X[i];
	}
	
	vector<Vertex2D<float32>> points;
	points = getEdgePoints<float32>(true, curve, 90, 0.08);

	//curve.shift_curve(true, p0, pn, tVec, gammaVec, angle1, angle2, alpha1m, alpha2);
	vector<float32> Pnx, Pny, Bnx, Bny, CurvY, CurvX, tx;
	for (size_t i = 0; i < 101; i++)
	{
		float32 t = float(i) / float(100);
		Bnx.push_back(curve.getPoint(t).x);
		Bny.push_back(curve.getPoint(t).y);
		CurvX.push_back(t);
		CurvY.push_back(curve.curvature(t));
	}
	for (auto &i : curve.PPoints)
	{
		Pnx.push_back(i.x);
		Pny.push_back(i.y);
	}
	//std::cout << fun.curvature_condition();
#ifndef _DEBUG
	//plt::figure();
	plt::subplot(2, 1, 1);
	plt::xlim(0, 1);
	plt::plot(Px, Py, "D--g");
	plt::plot(Bx, By);
	plt::plot(Pnx, Pny, "x--r");
	plt::plot(Bnx, Bny);
	vector<float32> x, y;
	x.push_back(point.x);
	y.push_back(point.y);
	plt::plot(x,y,"D");
	plt::plot(nearX, nearY, "X-r");
	vector<float32> px, py;
	px.push_back(points[0].x); px.push_back(points[1].x);
	py.push_back(points[0].y); py.push_back(points[1].y);
	plt::plot(px, py, "D");
	plt::grid(true);
	plt::axis("equal");
	plt::subplot(2, 1, 2);
	plt::grid(true);
	plt::xlim(0, 1);
	/*for (size_t i = 1; i < curve.PPoints.size() - 1; i++)
	{
		float32 t = curve.find_nearest(curve.PPoints[i]);
		tx.push_back(t);
	}*/
	//plt::plot(vector<float64>({ 0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0 }), f.curvature);
	plt::plot(CurvX, CurvY);
	plt::show();

#endif

    return 0;
	
}

