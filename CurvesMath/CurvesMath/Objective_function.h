#pragma once
#include "BezierCurve.h"


template <typename T>
class of_CamberLine
{
protected:
	BezierCurve<T>& assist_curve;
	vector<Vector2D<T>> tangent_vectors;
public:
	BezierCurve<T>& curve;
	Vertex2D<T>& point;
	T t;
	T betta, k1, k2;
	Vector2D<T> dt;
	vector<T> curvature;
	of_CamberLine() {};
	~of_CamberLine() {};
	of_CamberLine(BezierCurve<T> &_curve, BezierCurve<T> &_assist_curve) :curve(_curve), assist_curve(_assist_curve), point(Vertex2D<T>()), t(0.5), betta(0.0), dt(0.0, 0.0), k1(curve.curvature(0.0)), k2(curve.curvature(1.0)), curvature(10,100) {};
	of_CamberLine(BezierCurve<T> &_curve, BezierCurve<T> &_assist_curve, Vertex2D<T>& _point) :curve(_curve), assist_curve(_assist_curve), point(_point), t(0.5), betta(0.0), dt(0.0, 0.0), k1(curve.curvature(0.0)), k2(curve.curvature(1.0)), curvature(10, 100) { };
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
		k1 = curve.curvature(0.0);
		k2 = curve.curvature(1.0);
		curvatureRefresh(variables);
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
		k1 = curve.curvature(0.0);
		k2 = curve.curvature(1.0);
		curvatureRefresh(variables);
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
		for (size_t i = 2; i < size - 2; i++)
		{
			float32 t_assist = assist_curve.find_nearest(ppoints[i]);
			Vector2D<T> normal, tangent;
			tangent = assist_curve.dt(t_assist);
			normal.normal2vector(tangent, isUp);
			tempTangents.emplace_back(normal);
			Vertex2D<T> point_assist = assist_curve.getPoint(t_assist);
			T gamma_i = sqrtf(ppoints[i].length(point_assist));
			tempVariables.emplace_back(sqrtf(point_assist.x));
			tempVariables.emplace_back(gamma_i);
		}
		variables = tempVariables;
		tangent_vectors = tempTangents;
		recompute(variables);
		curvatureRefresh(variables);
	}
	T coincidence_condition(const vector<T>& variables)
	{
		Vertex2D<T> A = getBezierPoint(variables);
		return powf((A.x - point.x), 2) + powf((A.y - point.y), 2);
	}
	T tangent_condition()
	{
		return powf(dt.y - dt.x * tanf(betta*PI / 180.0), 2);
	}
	T curvature_condition()
	{
		float32 sum = 0.0;
		for (auto &i : curvature)
		{
			if (i <= 0)
				sum += 0;
			else
				sum += powf(EQUAL2EPS(i, 0.0, 0.8), 2);
		}
		return sum;
	}
	void curvatureRefresh(const vector<T> &variables)
	{
		curvature.clear();
		curvature.reserve(10);
		for (size_t i = 0; i < 11; i++)
		{
			float32 t = float(i) / 10.0;
			curvature.emplace_back(curve.curvature(t));
		}
	}
	T operator () (const vector<T>& variables)
	{
		curvatureRefresh(variables);
		auto B = coincidence_condition(variables) + tangent_condition() + curvature_condition();
		return (B);
	}
};

template<typename F, typename T>
class obj_function_alpha
{
public:

	obj_function_alpha() {};
	~obj_function_alpha() {};
	//obj_function_alpha(const of_CamberLine<T> &_f,const vector<T>& _var, const vector<T> &_p);
	//obj_function_alpha(of_CamberLine<T>& _f, vector<T>& _var, vector<T>& _p) :f(_f), var(_var), p(_p) {}
	obj_function_alpha(F &_f, vector<T> &_var, vector<T> &_p) : f(_f), var(_var), p(_p) {}
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

protected:
	const vector<T> &var, &p;
	F &f;
};

template<typename T>
class obj_functionCoorDescent
{
public:
	of_CamberLine<T> &f;
	int &index;
	obj_functionCoorDescent() {};
	~obj_functionCoorDescent() {};
	obj_functionCoorDescent(of_CamberLine<T>& _f, int &_index) :f(_f), index(_index) {};
	T operator () (vector<T> &variables)
	{
		return f(variables);
	}
	vector<Vertex2D<T>> getOriginalPPoints()
	{
		return f.curve.PPoints;
	}
private:
	//of_CamberLine<T> &f;
};

#ifndef _DEBUG

using namespace cppoptlib;
using Eigen::VectorXd;

template <typename T>
class CamberLineFunction : public Problem<double>
{
protected:
	BezierCurve<T>& assist_curve;
	vector<Vector2D<T>> tangent_vectors;
public:

	using typename cppoptlib::Problem<double>::Scalar;
	using typename cppoptlib::Problem<double>::TVector;

	BezierCurve<T>& curve;
	Vertex2D<T>& point;
	T t;
	T betta, k1, k2;
	Vector2D<T> dt;
	vector<T> curvature;
	CamberLineFunction() {};
	~CamberLineFunction() {};
	CamberLineFunction(BezierCurve<T> &_curve, BezierCurve<T> &_assist_curve)
		:curve(_curve), assist_curve(_assist_curve), point(Vertex2D<T>()), t(0.5),
		betta(0.0), dt(0.0, 0.0), k1(curve.curvature(0.0)), k2(curve.curvature(1.0)),
		curvature(11, 100) {};
	CamberLineFunction(BezierCurve<T> &_curve, BezierCurve<T> &_assist_curve, Vertex2D<T>& _point)
		:curve(_curve), assist_curve(_assist_curve), point(_point), t(0.5),
		betta(0.0), dt(0.0, 0.0), k1(curve.curvature(0.0)), k2(curve.curvature(1.0)),
		curvature(11, 100) {};
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
		k1 = curve.curvature(0.0);
		k2 = curve.curvature(1.0);
		curvatureRefresh(variables);
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
		k1 = curve.curvature(0.0);
		k2 = curve.curvature(1.0);

	}
	
	T coincidence_condition(const vector<T>& variables)
	{
		Vertex2D<T> A = getBezierPoint(variables);
		return powf((A.x - point.x), 2) + powf((A.y - point.y), 2);
	}
	T tangent_condition()
	{
		return powf(dt.y - dt.x * tanf(betta*PI / 180.0), 2);
	}
	T curvature_condition()
	{
		float32 sum = 0.0;
		for (auto &i : curvature)
		{
			if (i <= 0)
				sum += 0;
			else
				sum += powf(EQUAL2EPS(i, 0.0, 0.1), 2);
		}
		return sum;
	}
	void curvatureRefresh(const vector<T> &variables)
	{
		curvature.clear();
		curvature.reserve(11);
		for (size_t i = 0; i < 11; i++)
		{
			float32 t = float(i) / float(10);
			curvature.emplace_back(curve.curvature(t));
		}
	}
	double value(const TVector &x)
	{
		vector<T> variables;
		for (size_t i = 0; i < x.size(); i++)
		{
			variables.push_back(x[i]);
		}
		curvatureRefresh(variables);
		return coincidence_condition(variables) + tangent_condition() + curvature_condition();
	}
	void add_PPoint(TVector &x)
	{
		vector<T> tempVariables, variables;
		vector<Vector2D<T>> tempTangents;
		variables.reserve(x.size());
		variables.emplace_back(x[0]);
		variables.emplace_back(x[1]);
		tempVariables.reserve(variables.size() + 2);
		tempTangents.reserve(variables.size() - 2);
		curve.increase();
		vector<Vertex2D<float32>> &ppoints = curve.PPoints;
		size_t size = ppoints.size();
		bool isUp = true;
		if (curve.getPoint(0.5).y < assist_curve.getPoint(0.5).y)	isUp = false;	else isUp = true;

		tempVariables.emplace_back(sqrtf(ppoints[0].length(ppoints[1])));
		tempVariables.emplace_back(sqrtf(ppoints[size - 1].length(ppoints[size - 2])));
		for (size_t i = 2; i < size - 2; i++)
		{
			float32 t_assist = assist_curve.find_nearest(ppoints[i]);
			Vector2D<T> normal, tangent;
			tangent = assist_curve.dt(t_assist);
			normal.normal2vector(tangent, isUp);
			tempTangents.emplace_back(normal);
			Vertex2D<T> point_assist = assist_curve.getPoint(t_assist);
			T gamma_i = sqrtf(ppoints[i].length(point_assist));
			tempVariables.emplace_back(sqrtf(point_assist.x));
			tempVariables.emplace_back(gamma_i);
		}
		variables = tempVariables;
		x.resize(variables.size());
		for (size_t i = 0; i < variables.size(); i++)
		{
			x[i] = variables[i];
		}
		tangent_vectors = tempTangents;
		recompute(variables);
		curvatureRefresh(variables);

	}
	double operator = (const vector<T> &variables)
	{
		return coincidence_condition(variables) + tangent_condition() + curvature_condition();
	}
};
#endif // !_DEBUG




