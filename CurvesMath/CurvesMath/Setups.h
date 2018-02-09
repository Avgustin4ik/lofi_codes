#pragma once
#include <math.h>
#include <vector>
#include "Matrix.h"
#ifndef _DEBUG
#include "cppoptlib\meta.h"
#include "cppoptlib\problem.h"
#include "cppoptlib\solver\bfgssolver.h"
#include "cppoptlib\solver\lbfgsbsolver.h"
#include "cppoptlib\solver\newtondescentsolver.h"
#include "cppoptlib\solver\lbfgssolver.h"
#endif // !_DEBUG




typedef float float32;
typedef double float64;
const float32 INF = 0.0;
const float32 EPS = powf(10, -6);
const float32 MAX_SPLIT = 100;
typedef unsigned int uint;
const float64 PI = 3.141592653589793238462643;
const float32 LEFT_BORDER = 0.0;
const float32 RIGHT_BORDER = 1.0;

#define EQUAL2EPS(X,Y,DELTA) ABS(X-Y)<=DELTA ? Y : X

struct Configuration
{
public:
	const float32 curvature_eps;
	const uint iterations_limit;
	float32 h;//Приращение аргумента при нахождении производной	
	float32 alpha;//Коэффициент релаксации
	Configuration() :iterations_limit(200), h(0.001), alpha(0.1), curvature_eps(0.001) {};
	~Configuration() {};
};
