#pragma once
#include <math.h>
#include <vector>
#include "Matrix.h"


typedef float float32;
typedef double float64;
const float32 EPS = powf(10, -6);
const float32 MAX_SPLIT = 100;
typedef unsigned int uint;
const float64 PI = 3.141592653589793238462643;
const float32 LEFT_BORDER = 0.0;
const float32 RIGHT_BORDER = 1.0;



struct Configuration
{
public:
	const uint iterations_limit;
	float32 h;//Приращение аргумента при нахождении производной	
	float32 alpha;//Коэффициент релаксации
	Configuration():iterations_limit(50),h(0.001),alpha(0.1) {};
	~Configuration() {};
};
