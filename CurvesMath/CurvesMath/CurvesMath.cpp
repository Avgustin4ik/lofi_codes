// CurvesMath.cpp : Defines the entry point for the console application.
//
#include "stdafx.h"
#include "BezierCurve.h"
#include <iostream>
#include "MyMath.h"
#include "SidesFunction.h"
#include <time.h>
//************************************
#include "temp.h"
//************************************

#ifndef _DEBUG
#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;
using namespace cppoptlib;
using Eigen::VectorXd;
#endif

BezierCurve<float32> FishBone<float32>::skeleton;
int main()
{
	Configuration config;
	vector < Vertex2D<float32>>	P;
	P.push_back(Vertex2D<float32>(0, 0));
	P.push_back(Vertex2D<float32>(0.5 * cos(60 * PI / 180.0), 0.5 *sin(60 * PI / 180.0)));
	P.push_back(Vertex2D<float32>(1 - 0.5 *cos(21.5 * PI / 180.0), 0.5 *sin(21.5 * PI / 180.0)));
	P.push_back(Vertex2D<float32>(1, 0));
	BezierCurve<float32> curve(P);
	std::vector<float> x_(10);
	std::vector<float> y_(10);
	vector<float32> Bx(101), By(101);
	Vertex2D<float32> point(1, 0.5);
	Vertex2D<float32> nearPoint = curve.getPoint(curve.find_nearest(point));
	vector<float32> nearX, nearY;
	nearX.push_back(nearPoint.x);
	nearY.push_back(nearPoint.y);


	using POINTS = vector<Vertex2D<float32>>;
	using PP = Vertex2D<float32>;
	BezierCurve<float32> line_curve(POINTS({ PP(0.0,0.0), PP(1.0,0.0) }));
	//point = Vertex2D<float32>(0.16, 0.11);

	point = Vertex2D<float32>(0.25, 0.11);
	vector<float32> variables;
	variables.push_back(0.5);
	variables.push_back(0.5);
	of_CamberLine<float32> f(curve, line_curve, point);
	//minimization_Newthon<of_CamberLine<float32>, float32>(f, variables, config);
	//minimization_CoordinateDescent<of_CamberLine<float32>, float32>(F, variables, config);

#ifndef _DEBUG
	CamberLineFunction<float32> fun(curve, line_curve, point);
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
#endif // !_DEBUG
#ifdef _DEBUG
	curve.increase();
#endif	
	vector<Vertex2D<float32>> points_LEdge, points_TEdge;
	points_LEdge = getEdgePoints<float32>(true, curve, 120, 0.08);
	points_TEdge = getEdgePoints<float32>(false, curve, 90, 0.02);
	vector<float32> tVec, gammaVec;
	tVec.reserve(curve.PPoints.size() - 2);
	gammaVec.reserve(curve.PPoints.size() - 2);
	for (size_t i = 2; i < curve.PPoints.size() - 2; i++)
	{
		tVec.emplace_back(curve.find_nearest(curve.PPoints[i]));
		gammaVec.emplace_back(0.3);
	}
	float32 angle1 = 75;
	float32 angle2 = 30;
	BezierCurve<float32> suctionSide = curve.shift_curve(true, points_LEdge[0], points_TEdge[0], tVec, gammaVec, angle1, angle2, 0.5, 0.2);
	//	BezierCurve<float32> pressureSide = curve.shift_curve(false, points_LEdge[1], points_TEdge[1], tVec, gammaVec, angle1, angle2, 0.5, 0.2);

	float64 t_maxR = 0.18;
	float64 r_maxR = 0.1;
	Vertex2D<float32> centrMaxR = curve.getPoint(t_maxR);
	Circle<float32> maxR(r_maxR, centrMaxR);
	float64 t_bendR = 0.6;
	float64 r_bendR = 0.05;
	Vertex2D<float32> centrBend = curve.getPoint(t_bendR);
	Circle<float32> bend(r_bendR, centrBend);

	BoundaryConditions<float32> _conditions;
	VertexCondition<float32> _singeCondition;
	_conditions.data.reserve(4);
	_singeCondition.point = points_LEdge[0];
	_singeCondition.vector = Vector2D<float32>(75);
	_conditions.addCondition(_singeCondition);
	Vector2D<float32> _tv;
	_tv.normal2vector(curve.dt(t_maxR), true);
	_singeCondition.point = maxR.centr + _tv * maxR.radius;
	_singeCondition.vector = Vector2D<float32>().normal2vector(_tv, true);
	_conditions.addCondition(_singeCondition);
	_tv.normal2vector(curve.dt(t_bendR), true);
	_singeCondition.point = bend.centr + _tv * bend.radius;
	_singeCondition.vector = Vector2D<float32>().normal2vector(_tv, true);
	_conditions.addCondition(_singeCondition);
	_singeCondition.point = points_TEdge[0];
	_singeCondition.vector = Vector2D<float32>(-30);
	_conditions.addCondition(_singeCondition);
	variables.clear();
	variables.emplace_back(0.3);
	variables.emplace_back(0.3);
	variables.emplace_back(0.3);
	FishBone<float32>::skeleton = curve;
	SidesFunction<float32> functionSuctionSide(suctionSide, _conditions);
	BfgsSolver<SidesFunction<float32>> sSolver;
	//LbfgsSolver<SidesFunction<float32>> sSolver;
	//NewtonDescentSolver<SidesFunction<float32>> sSolver;
	stopCriteria.iterations = 100;
	stopCriteria.fDelta = 1e-3;
	stopCriteria.reset();
	sSolver.setStopCriteria(stopCriteria);
	VectorXd xx(3); xx << 0.5, 0.5, 0.5;
	//minimization_Newthon<SidesFunction<float32>, float32>(functionSuctionSide, variables, config);
#ifndef _DEBUG
	while (functionSuctionSide.value(xx) > 1e-3)
	{
		sSolver.minimize(functionSuctionSide, xx);
		cout << "f(x)		" << functionSuctionSide.value(xx) << "xx for SuctionSide		" << xx.transpose() << endl;
		vector<float32> v(xx.size());
		for (size_t i = 0; i < xx.size(); i++)	v[i] = xx[i];
		
		plt::clf();
		plt::subplot(3, 1, 1);
		plot_curve(functionSuctionSide.getCurve());
		plot_conditions(_conditions, 0.1);
		plot_fishbones(functionSuctionSide.getFishBones());
		plt::subplot(3, 1, 2);
		plot_vector2d(_conditions(0).vector);
		plot_vector2d(suctionSide.dt(suctionSide.find_nearest(_conditions(0).point)));
		plt::axis("equal");
		plt::subplot(3, 1, 3);
		plot_vector2d(_conditions(1).vector);
		plot_vector2d(suctionSide.dt(suctionSide.find_nearest(_conditions(1).point)));
		plt::axis("equal");
		plt::pause(0.1);
		if (functionSuctionSide.value(xx) > 1e-3)
		{
			functionSuctionSide.increaseCurve(v);
			xx.resize(v.size());
			for (size_t i = 0; i < xx.size(); i++)	xx[i] = v[i];
		}
	}

#endif // !_DEBUG
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
	vector < float32> Sx, Sy, PSx, PSy, Px, Py, PPx, PPy;
	for (size_t i = 0; i < 101; i++)
	{
		float32 t = float(i) / float(100);
		Sx.push_back(suctionSide.getPoint(t).x);
		Sy.push_back(suctionSide.getPoint(t).y);
		//		Px.push_back(pressureSide.getPoint(t).x);
		//		Py.push_back(pressureSide.getPoint(t).y);
	}
	for (auto &i : suctionSide.PPoints)
	{
		PSx.push_back(i.x);
		PSy.push_back(i.y);
	}
#ifndef _DEBUG
	plt::figure();
	plt::clf;
	plt::subplot(2, 2, 1);
	plt::xlim(0, 1);
	plt::plot(Pnx, Pny, "x--r");
	plt::plot(Bnx, Bny, "green");
	vector<float32> x, y;
	x.push_back(point.x);
	y.push_back(point.y);
	plt::plot(x, y, "D");
	plt::plot(nearX, nearY, "X-r");

	vector<float32> px, py, px2, py2;
	px.push_back(points_LEdge[0].x); px.push_back(points_LEdge[1].x);
	py.push_back(points_LEdge[0].y); py.push_back(points_LEdge[1].y);
	px2.push_back(points_TEdge[0].x); px2.push_back(points_TEdge[1].x);
	py2.push_back(points_TEdge[0].y); py2.push_back(points_TEdge[1].y);
	vector<float32> circleX(maxR.getXarray()), circleY(maxR.getYarray()), bendX(bend.getXarray()), bendY(bend.getYarray());

	plt::plot(circleX, circleY);
	plt::plot(bendX, bendY);
	plt::plot(px, py, "D");
	plt::plot(px2, py2, "D");
	plt::plot(Sx, Sy, "green");
	plt::plot(PSx, PSy, "X--b");

	plt::plot(Px, Py, "green");
	plt::plot(PPx, PPy, "X--b");

	plt::grid(true);
	plt::axis("equal");

	plt::subplot(2, 2, 2);
	plt::grid(true);
	plt::xlim(0, 1);
	plt::plot(CurvX, CurvY);

	//****График№3**********************
	plt::subplot(2, 2, 3);
	plt::axis("equal");
	plt::grid("true");
	plt::plot(Bnx, Bny, "green");
	x.push_back(point.x);
	y.push_back(point.y);
	plt::plot(x, y, "D");
	plt::plot(circleX, circleY);
	plt::plot(bendX, bendY);
	plt::plot(px, py, "D");
	plt::plot(px2, py2, "D");
	plt::plot(Sx, Sy, "green");
	plt::plot(PSx, PSy, "X--b");
	plt::plot(Px, Py, "green");
	plt::plot(PPx, PPy, "X--b");

	plot_conditions(_conditions, 0.1);

	//****График№3**********************

	plt::show();

#endif

	return 0;

}

