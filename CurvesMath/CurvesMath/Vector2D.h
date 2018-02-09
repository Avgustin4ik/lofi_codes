#pragma once
#include "Setups.h"
template <typename T>

struct Vector2D
{
public:
	T x, y;
	Vector2D() :x(0), y(0) {}
	~Vector2D() {}
	Vector2D(const T &_x, const T &_y) : x(_x), y(_y) { /*(*this).normalize(); */}
	Vector2D(const T &_angle) { x = cosf(_angle); y = sinf(_angle);}
//	Vector2D(T _x, T _y) : x(_x), y(_y) {}
	const float64 length() const
	{
		return sqrt(pow(x,2) + pow(y,2));
	}
	const Vector2D<T> normalize()
	{
		auto l = length();
		x /= l;
		y /= l;
		return *this;
	}
	const Vector2D<T> normal2vector(const Vector2D<T> _vec, bool isUp)
	{
		if (isUp) { x = -_vec.y; y = _vec.x; }
		else { x = _vec.y; y = -_vec.x; }
		return *this;
	}
	const T dot(const Vector2D<T> _vec)
	{
		return (x*_vec.x + y*_vec.y);
	}
	void reverse() 
	{
		this->x = this->x * -1;
		this->y = this->y * -1;
	}
	bool isHorizontal() { if (y = 0) return true };//Без погрешности
	bool isVertical() { if (x = 0) return true };//Без погрешности
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
	const Vector2D<T> operator + (const T &angle_rad)
	{
		return Vector2D(cos(atan(y / x) + angle_rad), sin(atan(y / x) + angle_rad));
	}
	const Vector2D<T> operator - (const T &angle_rad)
	{
		return Vector2D(cos(atan(y / x) - angle_rad), sin(atan(y / x) - angle_rad));
	}
private:

};
