#pragma once
#include "BezierCurve.h"
#include "Condition.h"
#ifndef _DEBUG

using namespace cppoptlib;
using Eigen::VectorXd;



template<typename T>
class SidesFunction : public Problem<double>
#else
template<typename T>
class SidesFunction
#endif
{
public:
#ifndef _DEBUG
	using typename cppoptlib::Problem<double>::Scalar;
	using typename cppoptlib::Problem<double>::TVector;
#endif // !_DEBUG

	SidesFunction();
	SidesFunction(BezierCurve<T> &_curve);
	SidesFunction(BezierCurve<T> &_curve, BoundaryConditions<T> &_conditions);
	SidesFunction(BezierCurve<T> &_curve, const vector<Vertex2D<T>> &_concidence_points, const vector<Vector2D<T>> &_tangent_vectors);
	~SidesFunction();
	void setConcidenceCondition(const vector<Vertex2D<T>> &_concidence_points);
	void setTangentCondition(const vector<Vector2D<T>> &_tangent_vectors);
	void increaseCurve(vector<T> &variables);
	void getVariables(vector<T> &variables);
	void setVariables(vector<T> &variables);
	BezierCurve<T>& getCurve();
	vector<Vertex2D<T>>& getConditionPoints();
	vector<FishBone<T>>& getFishBones();
	T operator () (vector<T> &variables);
#ifndef _DEBUG
	T operator () (const VectorXd &x);
	double value(const TVector &x);
#endif // _DEBUG

protected:

	vector<Vertex2D<T>> getPointsOnCurve();
	vector<Vector2D<T>> getTangentsOnCurve();
	void recompute_tangentsVector();
	void recompute_fishBones();
	T coincidenceCondition();
	T tangentCondition();
	T curvatureCondition();
	vector<T> tVec;					//Значение t для каждой точки, ближайшей к целевой
	vector<Vertex2D<T>> p;						//Целевые точки
	vector<Vector2D<T>> v;						//Вектора касательных в целевых точках
	vector<FishBone<T>> fishBones;
	vector<T> curvature;
	BezierCurve<T>& curve;
	//const BezierCurve<T> &camberLine;
	
};
template<typename T>
SidesFunction<T>::SidesFunction() 
	: curve(), tVec(), curvature(), p(), v()
{
}
template<typename T>
inline SidesFunction<T>::SidesFunction(BezierCurve<T>& _curve) :curve(_curve)//Незаконченный конструктор
{
}

template<typename T>
inline SidesFunction<T>::SidesFunction(BezierCurve<T> &_curve, BoundaryConditions<T> &_conditions)
	: tVec(), p(), v(), curvature(), fishBones(), curve(_curve)
{
	v.clear();
	v.reserve(_conditions.size()-2);
	p.reserve(_conditions.size()-2);
	fishBones.clear(); 
	fishBones.reserve(curve.PPoints.size() - 4);
	for (size_t i = 1; i < _conditions.size() - 1; i++)
	{
		p.emplace_back(_conditions(i).point);
		v.emplace_back(_conditions(i).vector);
		if ((i > 1) && (i < _conditions.size() - 1))	fishBones.emplace_back(FishBone<T>(curve.PPoints[i]));
	}
	tVec.clear(); tVec.reserve(p.size());
	for (auto &i : p)
		tVec.emplace_back(curve.find_nearest(i));
}

template<typename T>
inline SidesFunction<T>::SidesFunction(BezierCurve<T>& _curve, const vector<Vertex2D<T>>& _concidence_points, const vector<Vector2D<T>>& _tangent_vectors)
{
}

template<typename T>
inline void SidesFunction<T>::setConcidenceCondition(const vector<Vertex2D<T>> &_concidence_points) { p = _concidence_points; }

template<typename T>
inline void SidesFunction<T>::setTangentCondition(const vector<Vector2D<T>> &_tangent_vectors) { v = &_tangent_vectors; }

template<typename T>
inline void SidesFunction<T>::increaseCurve(vector<T>& variables)
{
	tVec.clear();
	variables.clear();
	curve.increase();
	vector<Vertex2D<float32>> &cP = curve.PPoints;
	variables.reserve(cP.size()-2);
	variables.emplace_back(sqrtf(cP[0].length(cP[1])));
	vector<Vertex2D<float32>>::iterator last = --(curve.PPoints.end());
	vector<Vertex2D<float32>>::iterator penultimate = last-1;
	variables.emplace_back(sqrtf(last->length(*penultimate)));
	recompute_fishBones();
	vector<FishBone<float32>>::iterator iterator_fb = fishBones.begin();
	for (size_t i = 2; i < cP.size() - 2; i++)
	{
		variables.emplace_back(iterator_fb->getValue());
		iterator_fb++;
	}
	recompute_tangentsVector();
}

template<typename T>
inline T SidesFunction<T>::operator()(vector<T>& variables)
{
	setVariables(variables);
	return coincidenceCondition() /*+ tangentCondition() + curvatureCondition()*/;
}
#ifndef _DEBUG
template<typename T>
inline T SidesFunction<T>::operator() (const VectorXd &x)
{
	vector<T> variables;
	for (size_t i = 0; i < x.size(); i++)	variables[i] = x[i];
	return this->operator()(variables);
}
template<typename T>
inline double SidesFunction<T>::value(const TVector & x)
{
	vector<T> variables;
	variables.reserve(x.size());
	for (size_t i = 0; i < x.size(); i++)
	{
		variables.emplace_back(x[i]);
	}
	setVariables(variables);
	return coincidenceCondition() + tangentCondition()/* + curvatureCondition()*/;
}
#endif // _DEBUG


template<typename T>
inline vector<Vertex2D<T>> SidesFunction<T>::getPointsOnCurve()
{

	vector<Vertex2D<T>> vertexVec;
	recompute_tangentsVector();
	vertexVec.reserve(p.size());
	for (auto &i : tVec)
		vertexVec.emplace_back(curve.getPoint(i));
	return vertexVec;
}

template<typename T>
inline vector<Vector2D<T>> SidesFunction<T>::getTangentsOnCurve()
{
	vector<Vector2D<T>> tangentsVec;
	tangentsVec.reserve(tVec.size());
	for (auto &i : tVec)
		tangentsVec.emplace_back(curve.dt(i));
	return tangentsVec;
}

template<typename T>
inline void SidesFunction<T>::getVariables(vector<T>& variables)
{
	vector<Vertex2D<T>> &cP = curve.PPoints;
	variables.clear();
	variables.resize(cP.size());
	variables.emplace_back(sqrtf(cP[0].length(cP[1])));
	variables.emplace_back(sqrtf(cP[cP.size() - 1].length(cP[cP.size() - 2])));
	for (size_t i = 2; i < cP.size() - 2; i++)
		variables.emplace_back(fishBones[i].getValue());
}

template<typename T>
inline void SidesFunction<T>::setVariables(vector<T>& variables)
{
	vector<Vertex2D<T>> &cP = curve.PPoints;
	Vector2D<T> inletTangent(curve.dt(0));
	Vector2D<T> outletTangent(curve.dt(1));
	if (outletTangent.x > 0)	outletTangent.reverse();
	cP[1] = cP[0] + inletTangent * powf(variables[0], 2);
	cP[cP.size() - 2] = cP[cP.size() - 1] + outletTangent * powf(variables[1], 2);
	for (size_t i = 2; i < cP.size() - 2; i++)
	{
		cP[i] = fishBones[i - 2].getPoint(variables[i]);
	} 
}

template<typename T>
inline BezierCurve<T>& SidesFunction<T>::getCurve()
{
	return curve;
}

template<typename T>
inline vector<Vertex2D<T>>& SidesFunction<T>::getConditionPoints()
{
	return p;
}

template<typename T>
inline vector<FishBone<T>>& SidesFunction<T>::getFishBones()
{
	return fishBones;
}

template<typename T>
inline void SidesFunction<T>::recompute_tangentsVector()
{
	tVec.clear();	tVec.reserve(p.size());
	for (auto &i : p)
		tVec.emplace_back(curve.find_nearest(i));
	tVec.shrink_to_fit();
}

template<typename T>
inline void SidesFunction<T>::recompute_fishBones()
{
	vector<Vertex2D<T>> &cP = curve.PPoints;
	fishBones.clear();	fishBones.reserve(cP.size() - 4);
	for (size_t i = 2; i < cP.size() - 2; i++)
	{
		fishBones.emplace_back(FishBone<T>(cP[i]));
	}
}

template<typename T>
inline T SidesFunction<T>::coincidenceCondition()
{
	vector<Vertex2D<T>> curvePointsVec(getPointsOnCurve());
	float64 coincidenceSum = 0;
	vector<Vertex2D<T>>::pointer ptr = &p[0];
	for (auto &i : curvePointsVec)
	{
		coincidenceSum += sqrtf(powf(i.x - (*ptr).x, 2) + powf(i.y - (*ptr).y, 2));
		ptr++;
	}
	return coincidenceSum;
}

template<typename T>
inline T SidesFunction<T>::curvatureCondition()
{
	return 0.0;
}
template<typename T>
inline T SidesFunction<T>::tangentCondition()
{
	vector<Vector2D<T>> tangents(getTangentsOnCurve());
	float64 resultSum = 0;
	for (size_t i = 0; i < tangents.size(); i++)
	{
		//resultSum += powf(tangents[i].x - v[i].x, 2) + powf(tangents[i].y - v[i].y, 2);
		resultSum += powf(tangents[i].y - tangents[i].x * tanf(v[i].y / v[i].x),2);
	}
	return resultSum;
}


template<typename T>
SidesFunction<T>::~SidesFunction()
{
}
