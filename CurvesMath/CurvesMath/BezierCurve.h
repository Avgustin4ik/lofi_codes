#pragma once
#include "Vertex2D.h"
#include "Vector2D.h"
#include <vector>
#include "Setups.h"


using namespace std;
template <typename T>

class BezierCurve
{
protected:
	
	float32 Bernstein(const uint m,const uint i)	//вычисление коэффициентов Ѕернштейна
	{
		float ans = float(fact(m)) / float(fact(m - i)) / float(fact(i));
	//	bernstein_data.push_back(ans);
		return ans;
	}
public:
	uint m;
//	vector<Vertex2D<T>> Points;	//координаты кривой
	vector<Vertex2D<T>> PPoints;	//координаты характерных точек
	vector<T> bernstein_data;

	BezierCurve()	//стандартна€ крива€ с 3 определ€ющими точками с координатами (0,0); (0.5,0.5); и (1,0)
		:m(2)
	{
		PPoints.push_back(Vertex2D<float32>(0, 0));
		PPoints.push_back(Vertex2D<float32>(0.5, 0.5));
		PPoints.push_back(Vertex2D<float32>(1, 0));
		compute();
	}
/*	BezierCurve(uint const _m)	//стандартна€ крива€ с 3 определ€ющими точками с координатами (0,0); (0.5,0.5); и (1,0)
		:m(_m)
	{
		size_t length = m + 1;
		for (size_t i = 0; i < length; i++)
		{
			PPoints.push_back(Vertex2D<T>(0, 0));
		}
		Compute();
	}*/
	BezierCurve(vector<Vertex2D<T>> _PPoints)
		:PPoints(_PPoints),m(_PPoints.size()-1)
	{
		compute();
	}
	BezierCurve(BezierCurve<T>& _curve) :PPoints(_curve.PPoints), m(_curve.m), bernstein_data(_curve.bernstein_data) {};
	float32 fact(const float value) {	//вычисление факториала числа типа float
		float ans = value;
		if (value == 0) { ans = 1; return ans; }

		if (value < 0)
			return 0;
		if (value > 0)
		{
			for (int i = 2; i < value; i++)
				ans *= i;
			return ans;
		}
	}
	void compute() {
		bernstein_data.clear();
			for (size_t i = 0; i < PPoints.size(); i++)
			{
					auto b = Bernstein(m, i);
//					R.x += B*PPoints[i].x;
//					R.y += B*PPoints[i].y;
					bernstein_data.push_back(b);
			}
		//	Points.push_back(R);
	}
	const Vertex2D<T> getPoint(const float32 t)
	{
		if ((t < 0) || (t > 1))	throw(exception());
		Vertex2D<T> R(0, 0);
		for (size_t i = 0; i < PPoints.size(); i++)
		{
			R = R + PPoints[i] * bernstein_data[i] * (powf(t, i)*powf(1 - t, m - i));
//			R.x += bernstein_data[i] * (powf(t, i)*powf(1 - t, m - i))*PPoints[i].x;
//			R.y += bernstein_data[i] * (powf(t, i)*powf(1 - t, m - i)) * PPoints[i].y;
//			R.x += Bernstein(m, i, t)*PPoints[i].x;
//			R.y += Bernstein(m, i, t)*PPoints[i].y;
		}
		return R;
	}
	const Vertex2D<T> operator [] (const float32 t)
	{
		if ((t < 0) || (t > 1))	throw(exception());
		Vertex2D<T> R(0, 0);
		for (size_t i = 0; i < PPoints.size(); i++)
		{
			R.x += bernstein_data[i] * (powf(t, i)*powf(1 - t, m - i)) * PPoints[i].x;
			R.y += bernstein_data[i] * (powf(t, i)*powf(1 - t, m - i)) * PPoints[i].y;
		//	R.x += Bernstein(m, i, t)*PPoints[i].x;
		//	R.y += Bernstein(m, i, t)*PPoints[i].y;
		}
		return R;
	}
	const Vector2D<T> dt(const float32 t)
	{
		if ((t < 0) || (t > 1))	throw(exception());
		Vector2D<T> ans;
		for (size_t i = 0; i < PPoints.size() - 1; i++)
		{
			ans.x += m*Bernstein(m - 1, i)*(powf(t, i)*powf(1 - t, m - 1 - i))*(PPoints[i + 1].x - PPoints[i].x);
			ans.y += m*Bernstein(m - 1, i)*(powf(t, i)*powf(1 - t, m - 1 - i))*(PPoints[i + 1].y - PPoints[i].y);
		}
		return ans.normalize();
	}
	const Vector2D<T> ddt(const float32 t)
	{
		if ((t < 0) || (t > 1))	throw(exception());
		Vector2D<T> ans;
		for (size_t i = 0; i < PPoints.size() - 2; i++)
		{
			ans.x += m*(m - 1)*Bernstein(m - 2, i)*(powf(t, i)*powf(1 - t, m - 2 - i))*(PPoints[i + 2].x - 2 * PPoints[i + 1].x + PPoints[i].x);
			ans.y += m*(m - 1)*Bernstein(m - 2, i)*(powf(t, i)*powf(1 - t, m - 2 - i))*(PPoints[i + 2].y - 2 * PPoints[i + 1].y + PPoints[i].y);
		}
		return ans.normalize();
	}
	void increase()
	{
		vector<Vertex2D<T>> newPPoints;
		newPPoints.push_back(PPoints[0]);
		m++;
		for (size_t i = 1; i < m; i++)
		{
			Vertex2D<T> R;
			R.x = float(i) / float(m)*PPoints[i - 1].x + (float(1) - float(i) / float(m))*PPoints[i].x;
			R.y = float(i) / float(m)*PPoints[i - 1].y + (float(1) - float(i) / float(m))*PPoints[i].y;
			newPPoints.push_back(R);
		}
		newPPoints.push_back(PPoints.back());
		PPoints.swap(newPPoints);
		compute();
	}
	const T find_nearest(const Vertex2D<T>& point)
	{
		function_bezier2point<float32> f(*this,point);
		vector<T> variables;
		T t = 0.5;
		variables.push_back(t);
		method_bisection(f, variables, LEFT_BORDER, RIGHT_BORDER);
//		newton_minimization_simple(f, variables, LEFT_BORDER, RIGHT_BORDER);
		t = variables.at(0);
		return (t);
	}
	const Vertex2D<T> operator () (const float32 t) const
	{
		if ((t < 0) || (t > 1))	throw(exception());
		Vertex2D<T> R(0, 0);
		for (size_t i = 0; i < PPoints.size(); i++)
		{
			R.x += bernstein_data[i] * (powf(t, i)*powf(1 - t, m - i)) * PPoints[i].x;
			R.y += bernstein_data[i] * (powf(t, i)*powf(1 - t, m - i)) * PPoints[i].y;
			//	R.x += Bernstein(m, i, t)*PPoints[i].x;
			//	R.y += Bernstein(m, i, t)*PPoints[i].y;
		}
		return R;
	}
	const T curvature(const T t)
	{
		auto ans = (dt(t).x*ddt(t).y - dt(t).y*ddt(t).x) / powf((powf(dt(t).x, float(2)) + powf(dt(t).y, float(2))), float(3) / float(2));
		return ans;
	}
	BezierCurve<T> shift_curve(const bool isSuctionSide, const Vertex2D<T> _p0, const Vertex2D<T> _pn, const vector<T> _t_arr, const vector<T> _gamma, T _angle1, T _angle2, T _alpha1, T _alpha2)
	{
		vector<Vertex2D<T>> newPoints;
		newPoints.push_back(_p0);
		Vector2D<float32> vtemp_1(_angle1);
		Vector2D<float32> vtemp_n(_angle2);
		Vertex2D<float32> p1;
		p1 = newPoints[0] + vtemp_1*_alpha1;
		Vertex2D<float32> pn_1 = _pn + vtemp_n*_alpha2;
		newPoints.push_back(p1);
		for (size_t i = 0; i < _t_arr.size(); i++)
		{
			Vertex2D<float32> ptemp;
			ptemp =  getPoint(_t_arr[i]);
			Vector2D<float32> vtemp;
			vtemp.normal2vector(dt(_t_arr[i]),isSuctionSide);
			vtemp.normalize();
			vtemp = vtemp * _gamma[i];
			ptemp = ptemp + vtemp;
			newPoints.push_back(ptemp);
		}
		newPoints.push_back(pn_1);
		newPoints.push_back(_pn);
		return BezierCurve<T>(newPoints);
	}
	
};
