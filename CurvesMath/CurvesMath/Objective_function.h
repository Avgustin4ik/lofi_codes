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
class objective_function_tangent
{
protected:
	BezierCurve<T>& assist_curve;
	vector<Vector2D<T>> tangent_vectors;
public:
	BezierCurve<T>& curve;
	Vertex2D<T>& point;
	T t;
	T betta;
	Vector2D<T> dt;
	objective_function_tangent() {};
	~objective_function_tangent() {};
	objective_function_tangent(BezierCurve<T> &_curve, BezierCurve<T> &_assist_curve) :curve(_curve), assist_curve(_assist_curve), point(Vertex2D<T>()), t(0.5), betta(0.0), dt(0.0, 0.0) {};
	objective_function_tangent(BezierCurve<T> &_curve, BezierCurve<T> &_assist_curve, Vertex2D<T>& _point) :curve(_curve),assist_curve(_assist_curve), point(_point), t(0.5), betta(0.0), dt(0.0,0.0) { };
	Vertex2D<T> getBezierPoint(const vector<T>& variables)
	{
		vector<Vertex2D<T>> &ppoints = curve.PPoints;
		vector<Vertex2D<T>> temp_ppoints = curve.PPoints;
		T x1(variables[0]);
		T x2(variables[1]);
		Vector2D<T> temp = curve.dt(0);
		if (temp.y < 0 | temp.x < 0)	temp.reverse();
		size_t last = ppoints.size() - 1;
		Vertex2D<T>& P1 = ppoints[1];
		Vertex2D<T>& P2 = ppoints[last - 1];
		P1 = ppoints[0] + temp*powf(x1, 2);
		temp = curve.dt(1);//вот тут вопросы!!!!!!!
		if (temp.x > 0 | temp.y < 0)	temp.reverse();
		P2 = ppoints[last] + temp*powf(x2, 2);
		if (variables.size() > 2)
		{
			size_t j = 2;
			for (size_t i = 2; i < variables.size() - 1; i += 2)
			{
				ppoints[j] = Vertex2D<T>(powf(variables[i], 2), powf(variables[i + 1], 2));
				j++;
			}
		}
		auto tt = curve.find_nearest(point);
		dt = curve.dt(tt);
		auto R = curve.getPoint(tt);
		ppoints = temp_ppoints;
		return R;
	}

	void recompute(const vector<T>& variables)
	{
		vector<Vertex2D<T>> &ppoints = curve.PPoints;
		T x1(variables[0]);
		T x2(variables[1]);
		Vector2D<T> temp = curve.dt(0);
		if (temp.y < 0 | temp.x < 0)	temp.reverse();
		size_t last = ppoints.size() - 1;
		Vertex2D<T>& P1 = ppoints[1];
		Vertex2D<T>& P2 = ppoints[last - 1];
		P1 = ppoints[0] + temp*powf(x1, 2);
		temp = curve.dt(1);//вот тут вопросы!!!!!!!
		if (temp.x > 0 | temp.y < 0)	temp.reverse();
		P2 = ppoints[last] + temp*powf(x2, 2);
		if (variables.size() > 2)
		{
			size_t j = 2;
			for (size_t i = 2; i < variables.size() - 1; i += 2)
			{
				ppoints[j] = Vertex2D<T>(powf(variables[i], 2), powf(variables[i + 1], 2));
				j++;
			}
		}
			
		t = curve.find_nearest(point);
		dt = curve.dt(t);
	}
	void add_PPoint(vector<T>& variables)
	{
		vector<T> tempVariables;
		vector<Vector2D<T>> tempTangents;
		tempVariables.reserve(variables.size() + 2);
		tempTangents.reserve(variables.size() - 2);
		curve.increase();
		vector<Vertex2D<float32>> &ppoints = curve.PPoints;
		size_t size = ppoints.size();
		bool isUp = true;
		if (curve.getPoint(0.5).y < assist_curve.getPoint(0.5).y)	isUp = false;	else isUp = true;

		tempVariables.emplace_back(sqrtf(ppoints[0].length(ppoints[1])));
		tempVariables.emplace_back(sqrtf(ppoints[size - 1].length(ppoints[size - 2])));
		for (size_t i = 2; i < size-2; i++)
		{
			float32 t_assist = assist_curve.find_nearest(ppoints[i]);
			Vector2D<T> normal, tangent;
			tangent = assist_curve.dt(t_assist);
			normal.normal2vector(tangent,isUp);
			tempTangents.emplace_back(normal);
			Vertex2D<T> point_assist = assist_curve.getPoint(t_assist);
			T gamma_i = sqrtf(ppoints[i].length(point_assist));
			tempVariables.emplace_back(sqrtf(point_assist.x));
			tempVariables.emplace_back(gamma_i);
		}
		variables = tempVariables;
		tangent_vectors = tempTangents;
		recompute(variables);
	}
	T coincidence_condition()
	{
		
	}
	T tangent_condition()
	{

	}
	T operator () (const vector<T>& variables)
	{
		Vertex2D<T> A = getBezierPoint(variables);
		auto B = powf((A.x - point.x), 2) + powf((A.y - point.y), 2) + powf(dt.y - dt.x * tanf(betta*PI / 180.0), 2);
		return (B);
	}

};

template<typename T>
class obj_function_alpha
{
public:
	obj_function_alpha() {};
	~obj_function_alpha() {};
	//obj_function_alpha(const objective_function_tangent<T> &_f,const vector<T>& _var, const vector<T> &_p);
	obj_function_alpha(objective_function_tangent<T>& _f, vector<T>& _var, vector<T>& _p) :f(_f), var(_var), p(_p) {}
	T operator()(vector<T>& alpha)
	{
		vector<T> variables;
		variables.reserve(var.size());
		for (size_t i = 0; i < var.size(); i++)
		{
			variables.emplace_back(var[i] + alpha[0] * p[i]);
		}

		return f(variables);
	}
private:
	const vector<T> &var, &p;
	objective_function_tangent<T> &f;
};

template<typename T>
class obj_functionCoorDescent
{
public:
	int &index;
	obj_functionCoorDescent() {};
	~obj_functionCoorDescent() {};
	obj_functionCoorDescent(objective_function_tangent<T>& _f, int &_index) :f(_f), index(_index) {};
	T operator () (vector<T> &variables)
	{
		return f(variables);
	}
	vector<Vertex2D<T>> getOriginalPPoints()
	{
		return f.curve.PPoints;
	}
private:
	objective_function_tangent<T> &f;
};

