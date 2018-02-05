#pragma once
#include "Objective_function.h"
template<typename T>

class Curvature:public of_CamberLine<T>
{
public:
	Curvature();
	~Curvature();

	T operator () (const T _t);
private:

};

template<typename T>
inline Curvature<T>::Curvature()
{
}

template<typename T>
inline Curvature<T>::~Curvature()
{
}

template<typename T>
inline T Curvature<T>::operator()(const T _t)
{
	
	return T();
}
