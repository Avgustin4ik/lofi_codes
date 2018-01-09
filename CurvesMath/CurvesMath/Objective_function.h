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
	objective_function(BezierCurve<T> &_curve) :curve(_curve), point(Vertex2D<T>()) {};
	objective_function(BezierCurve<T> &_curve, Vertex2D<T>& _point) :curve(_curve), point(_point), t(0.5) { };
	Vertex2D<T> getBezierPoint(const T& x1, const T& x2) const
	{
		if ((t < 0) || (t > 1))	throw(exception());
		Vertex2D<T> R(0, 0);
		Vector2D<T> temp = curve.dt(0);
		if (temp.y < 0 | temp.x < 0)	temp.reverse();
		Vertex2D<T>& P1 = curve.PPoints[1];
		Vertex2D<T>& P2 = curve.PPoints[2];
		P1 = curve.PPoints[0] + temp*x1;
		temp = curve.dt(1);//вот тут вопросы!!!!!!!
		if (temp.x > 0 | temp.y < 0)	temp.reverse();
		P2 = curve.PPoints[3] + temp*x2;
		auto tn = curve.find_nearest(point);
		for (size_t i = 0; i < curve.PPoints.size(); i++)
		{
			R = R + curve.PPoints[i] * curve.bernstein_data[i] * (powf(tn, i)*powf(1 - tn, curve.m - i));
		}
		return R;
	}
	T operator () (const vector<T>& var_arr) const
	{
		if ((t < 0) || (t > 1))	throw(exception());
		T x1(var_arr[0]), x2(var_arr[1]);
		Vertex2D<T> A = getBezierPoint(x1, x2);
		auto B = A.length(point);
		return (B);
	}

};
