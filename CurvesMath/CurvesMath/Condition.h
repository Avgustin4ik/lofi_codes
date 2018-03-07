#pragma once
#include "Setups.h"
#include "Vertex2D.h"

template<typename T>
struct VertexCondition
{
public:
	Vertex2D<T> point;
	Vector2D<T> vector;
};
template<typename T>
struct BoundaryConditions
{
public:
	std::vector<VertexCondition<T>> data;
	size_t size()
	{
		return data.size();
	}
	void addCondition(const Vertex2D<T> _point, const Vector2D<T> _vector)
	{
		VertexCondition<T> c;
		c.point(_point);
		c.vector(_vector);
		data.emplace_back(c);
	}
	void addCondition(const VertexCondition<T> _vc)
	{
		data.emplace_back(_vc);
	}
	VertexCondition<T> operator ()(const size_t &i)
	{
		return data[i];
	}
};