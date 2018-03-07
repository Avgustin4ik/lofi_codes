#pragma once
#include "Setups.h"
#include "MyMath.h"

template<typename T>
class FishBone
{
public:
	//static BezierCurve<T> &skeleton;
	static BezierCurve<T> skeleton;
	FishBone();
	FishBone(const float64 &tSk, const T &variable);
	FishBone(const Vertex2D<T> &endPoint);
	~FishBone();

	const BezierCurve<T>& getSkeleton();
	const Vertex2D<T> getPoint(const T& x);
	const Vertex2D<T> getPoint();
	const Vertex2D<T> getPointOnSkeleton();
	void operator()(const T &x);
	void recompute(Vertex2D<T> &p, T &gamma);
	const T getValue();

	bool isSuctionSide;
protected:
	Vertex2D<T> skeletonPoint, endPoint;
	Vector2D<T> bone;
	
};
template<typename T>
FishBone<T>::FishBone()
	:skeletonPoint(), endPoint(), bone(), isSuctionSide(true)
{
}
template<typename T>
inline FishBone<T>::FishBone(const float64 &tSk, const T &variable)
	:skeletonPoint(skeleton.getPoint(tSk)), bone(skeleton.dt(tSk)), isSuctionSide(true)
{
	bone.normalize();
	endPoints = skeletonPoint + bone * powf(var, 2);
}
template<typename T>
inline FishBone<T>::FishBone(const Vertex2D<T> &_endPoint)
	:endPoint(_endPoint), isSuctionSide(true)
{
	auto t = skeleton.find_nearest(endPoint);
	Vector2D<T> dt = skeleton.dt(t);
	bone.normal2vector(dt, isSuctionSide);
	bone.normalize();
	skeletonPoint = skeleton.getPoint(t);
}
template<typename T>
FishBone<T>::~FishBone()
{
}

template<typename T>
inline const BezierCurve<T>& FishBone<T>::getSkeleton()
{
	return &skeleton;
}

template<typename T>
inline const Vertex2D<T> FishBone<T>::getPoint(const T & x)
{
	endPoint = skeletonPoint + bone * powf(x, 2);
	return endPoint;
}

template<typename T>
inline const Vertex2D<T> FishBone<T>::getPoint()
{
	
	// TODO: вставьте здесь оператор return
	return endPoint;
}

template<typename T>
inline const Vertex2D<T> FishBone<T>::getPointOnSkeleton()
{
	return skeletonPoint;
}

template<typename T>
inline void FishBone<T>::operator()(const T & x)
{
	float64 t = skeleton.find_nearest(skeletonPoint);
	bone = skeleton.dt(t);
	endPoint = skeletonPoint + bone * powf(x, 2);
}

template<typename T>
inline void FishBone<T>::recompute(Vertex2D<T> &p, T &gamma)
{
	float64 t = skeleton.find_nearest(p);
	bone = skeleton.dt(t).normalize();
	skeletonPoint = skeleton.getPoint(t);
	endPoint = skeletonPoint + bone * powf(gamma, 2);
}

template<typename T>
inline const T FishBone<T>::getValue()
{
	return sqrtf(skeletonPoint.length(endPoint));
}