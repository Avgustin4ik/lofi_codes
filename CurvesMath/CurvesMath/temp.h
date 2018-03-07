#pragma once
#include "SidesFunction.h"
#include "Condition.h"
#ifndef _DEBUG
#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;
void plot_conditions(const BoundaryConditions<float32> &_conditions, const float32 &_scale)
{
	vector<double> x, y, vx, vy;
	for (auto &i : _conditions.data)
	{
		x.push_back(i.point.x);
		y.push_back(i.point.y);
		vx.push_back(i.point.x);
		vy.push_back(i.point.y);
		Vector2D<float32> temp = i.vector;
		Vertex2D<float32> end(i.point + temp * _scale);
		vx.push_back(end.x);
		vy.push_back(end.y);
		plt::plot(vx, vy, "-r");
		vx.clear();
		vy.clear();
	}
	plt::plot(x, y, "o");
}

void plot_curve(BezierCurve<float32> &_curve)
{
	vector<float64> Bx, By, Px, Py;
	for (auto i = 0; i < 101; i++) {
		float64 z = float(i) / float(100);
		Bx.push_back(_curve.getPoint(z).x);
		By.push_back(_curve.getPoint(z).y);
	}
	for (auto &i : _curve.PPoints) {
		Px.push_back(i.x);
		Py.push_back(i.y);
	}
	plt::grid(true);
	plt::axis("equal");
	plt::plot(Px, Py, "x--r");
	plt::plot(Bx, By);
}
void plot_fishbones(vector<FishBone<float32>> &bones)
{
	for (auto &i : bones)
	{
		vector<float64> x, y;
		Vertex2D<float32> a(i.getPoint());
		Vertex2D<float32> b(i.getPointOnSkeleton());
		x.push_back(a.x);	x.push_back(b.x);
		y.push_back(a.y);	y.push_back(b.y);
		plt::plot(x, y, "s-.y");
	}
}

void plot_vector2d(const Vector2D<float32> &v)
{
	vector<float64> x, y;
	x.push_back(0.0);	x.push_back(v.x);
	y.push_back(0.0);	y.push_back(v.y);
	plt::plot(x, y);
}
#endif // !_DEBUG

