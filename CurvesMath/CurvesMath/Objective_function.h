#pragma once
#include "BezierCurve.h"
template <typename T>

class objective_function
{
public:
	BezierCurve<T>& curve;
	Vertex2D<T>& point;
	T t;
	objective_function() {};
	~objective_function() {};
	objective_function(BezierCurve<T> &_curve) :curve(_curve), point(Vertex2D<T>()),t(0.5) {};
	objective_function(BezierCurve<T> &_curve, Vertex2D<T>& _point) :curve(_curve), point(_point), t(0.5) { };
	Vertex2D<T> getBezierPoint(const T& x1, const T& x2)
	{
		Vector2D<T> temp = curve.dt(0);
		if (temp.y < 0 | temp.x < 0)	temp.reverse();
		auto p1 = curve.PPoints[1];
		auto p2 = curve.PPoints[2];
		Vertex2D<T>& P1 = curve.PPoints[1];
		Vertex2D<T>& P2 = curve.PPoints[2];
		P1 = curve.PPoints[0] + temp*x1;
		temp = curve.dt(1);//вот тут вопросы!!!!!!!
		if (temp.x > 0 | temp.y < 0)	temp.reverse();
		P2 = curve.PPoints[3] + temp*x2;
		T tt = curve.find_nearest(point);
		if ((tt < 0) || (tt > 1))	throw(exception());
		auto R = curve.getPoint(tt);
		P1 = p1;
		P2 = p2;
		return R;
	}
	void recompute(const T& x1, const T& x2)
	{
		Vector2D<T> temp = curve.dt(0);
		if (temp.y < 0 | temp.x < 0)	temp.reverse();
		auto p1 = curve.PPoints[1];
		auto p2 = curve.PPoints[2];
		Vertex2D<T>& P1 = curve.PPoints[1];
		Vertex2D<T>& P2 = curve.PPoints[2];
		P1 = curve.PPoints[0] + temp*x1;
		temp = curve.dt(1);//вот тут вопросы!!!!!!!
		if (temp.x > 0 | temp.y < 0)	temp.reverse();
		P2 = curve.PPoints[3] + temp*x2;
		t = curve.find_nearest(point);
	}
	T operator () (const vector<T>& var_arr)
	{
		T x1(var_arr[0]), x2(var_arr[1]);
		Vertex2D<T> A = getBezierPoint(x1,x2);
		auto B = powf((A.x - point.x), 2) + powf((A.y - point.y), 2);
		return (B);
	}

};

template <typename T>
class obj_function_tangent_point
{
public:
	BezierCurve<T>& curve;
	Vertex2D<T>& point;
	T t;
	obj_function_tangent_point() {};
	~obj_function_tangent_point() {};
	obj_function_tangent_point(BezierCurve<T> &_curve) :curve(_curve), point(Vertex2D<T>()), t(0.5) {};
	obj_function_tangent_point(BezierCurve<T> &_curve, Vertex2D<T>& _point) :curve(_curve), point(_point), t(0.5) { };
	void recompute(const T& x1, const T& x2)
	{
		Vector2D<T> temp = curve.dt(0);
		if (temp.y < 0 | temp.x < 0)	temp.reverse();
		auto p1 = curve.PPoints[1];
		auto p2 = curve.PPoints[2];
		Vertex2D<T>& P1 = curve.PPoints[1];
		Vertex2D<T>& P2 = curve.PPoints[2];
		P1 = curve.PPoints[0] + temp*x1;
		temp = curve.dt(1);//вот тут вопросы!!!!!!!
		if (temp.x > 0 | temp.y < 0)	temp.reverse();
		P2 = curve.PPoints[3] + temp*x2;
		t = curve.find_nearest(point);
	}
	T operator () (const vector<T>& var_arr)
	{
		T x1(var_arr[0]), x2(var_arr[1]);
		recompute(x1, x2);
		return (curve.dt(t).x);
	}


};