#pragma once
#include "Test.h"
template <typename T>
struct Vertex2D
{
public:
	T x, y;
	Vertex2D() :x(0), y(0) {}
	Vertex2D(T _x, T _y) :x(_x), y(_y) {}
	void operator = (Vertex2D<T> point)
	{
		x = point.x;
		y = point.y;
	}
	Vertex2D<T> operator-(const Vertex2D<T> &point) const
	{
		Vertex2D<T> R;
		R.x = x - point.x;
		R.y = y - point.y;
		return R;
	}
	Vertex2D<T> operator/(const T &val) const
	{
		Vertex2D<T> R;
		R.x = x/val;
		R.y = y/val;
		return R;
	}
	const T length(const Vertex2D<T> &point)
	{
		auto R = sqrtf(powf((x - point.x), 2) + powf((y - point.y), 2));
		return R;
	}
};


