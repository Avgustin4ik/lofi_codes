#pragma once
#include "Test.h"
template <typename T>

struct Vector2D
{
public:
	T x, y;
	Vector2D() :x(0), y(0) {}
	~Vector2D() {}
	Vector2D(const T &_x, const T &_y) : x(_x), y(_y) {}
//	Vector2D(T _x, T _y) : x(_x), y(_y) {}
	const float64 length() const
	{
		return sqrt(pow(x,2) + pow(y,2));
	}
	const void normalize()
	{
		auto l = length();
		x /= l;
		y /= l;
	}
	const T operator dot(const Vector2D<T> _vec)
	{
		return (x*_vec.x + y*_vec.y);
	}
	Vector2D<T> operator = (const Vector2D<T> _v)
	{
		x = _v.x;
		y = _v.y;
		return *this;
	}
	Vector2D<T> operator *(const T &_value)
	{
		T X = x*_value;
		T Y = y*_value;
		return Vector2D(X, Y);
	}
	Vector2D<T> operator *=(const T &_value)
	{
		T X = x*_value;
		T Y = y*_value;
		return Vector2D(X, Y);
	}
	const Vector2D<T> operator +(const Vector2D<T> _vec)
	{
		T X = x + _vec.x;
		T Y = y + _vec.y;
		return Vector2D(X, Y);
	}
private:

};
