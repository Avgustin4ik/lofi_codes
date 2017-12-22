// CurvesMath.cpp : Defines the entry point for the console application.
//
#include "stdafx.h"
#include "BezierCurve.h"
#include <iostream>
#include "MyMath.h"
#include <time.h>



int main()
{
	Vertex2D<float32> a(0,0);
	Vertex2D<float32> b(1, 1);
	Vertex2D<float32> c(2, 0);
	vector<Vertex2D<float32>> ppoints;
	ppoints.push_back(a);
	ppoints.push_back(b);
	ppoints.push_back(c);
	BezierCurve<float32> curve(ppoints);
	auto P0 = Vertex2D<float32>(-1, 1);
	auto PN = Vertex2D<float32>( 3, 1);
	vector<float32> T_value;
	T_value.push_back(0.33);
	T_value.push_back(0.5);
	T_value.push_back(0.66);
	vector<float32> Gamma_value;
	Gamma_value.push_back(2);
	Gamma_value.push_back(4);
	Gamma_value.push_back(3);
	float32 Angle1 = 45;
	float32 Angle2 = 90+45;
	float32 Alpha1 = 1;
	float32 Alpha2 = 1;
	auto A = curve.getPoint(0.5);
	bool isSuctionSide = true;
	BezierCurve<float32> curve2 = curve.shift_curve(isSuctionSide,P0, PN, T_value, Gamma_value, Angle1, Angle2, Alpha1, Alpha2);
	BezierCurve<float32> curve3 = curve.shift_curve(bool(false), Vertex2D<float32>(1,-1), Vertex2D<float32>(3, -1), T_value, Gamma_value, Angle1, Angle2, Alpha1, Alpha2);
//	function<float32> F;
/*	float32 x1 = 0, x2 = 0, x1n = 1, x2n = 1;
	float32 h = 0.001;
	using d_f = derivarive<function<float32>, float32>;
	using dd_f = second_derivative<function<float32>, float32>;
	d_f dF(F,h);
	dd_f ddF(F, h);
	vector<float32> var_arr;
	var_arr.push_back(x1n);
	var_arr.push_back(x2n);
	Vertex2D<float32> VX(x1, x2);
	Vertex2D<float32> VXN(x1n, x2n);
	uint iterator = 0;
	while (sqrtf(powf(x1-x1n,2) + powf(x2-x2n,2)) > EPS)
	{
		iterator++;
		x1 = x1n;
		x2 = x2n;
		
		auto g1 = dF(var_arr, 0,t);
		auto g2 = dF(var_arr, 1);
		auto G1 = ddF(var_arr, 0, 0);
		auto G2 = ddF(var_arr, 0, 1);
		auto G3 = ddF(var_arr, 1, 0);
		auto G4 = ddF(var_arr, 1, 1);
		float64 p1, p2;

		float64 values[2][3] = {
			{ G1, G3, -g1 },
			{ G2, G4, -g2 }
		};
		slau(values, p1, p2);
		x1n += p1;
		x2n += p2;
		VXN.x = x1n;
		VXN.y = x2n;
		var_arr[0] = x1n;
		var_arr[1] = x2n;
	}*/
	float32 x1 = 0.5, x2 = 0.5;
	float32 h = 0.001;
	vector<float32> var_arr;
	var_arr.push_back(x1);
	var_arr.push_back(x2);
	curve.increase();
	Vertex2D <float32> point(0.5,0.5);
	objective_function<float32> obj_function(curve,point);
	obj_function.t = 0.5;
	auto AA  =  obj_function.getBezierPoint(0.5, 1.5);
//	curve.find_nearest(point);
//	auto D = newton_minimization(obj_function, var_arr);
//	*******************************************
	temp_function<float32> f;
	float32 x = 0.5;
	var_arr.clear();
	var_arr.push_back(x);
	newton_minimization_simple(f, var_arr, LEFT_BORDER, RIGHT_BORDER);
	c.y = 1;
	Vertex2D<float32> d(3, 0);
	vector < Vertex2D<float32>>	P;
	P.push_back(a);
	P.push_back(b);
	P.push_back(c);
	P.push_back(d);
	BezierCurve<float32> CURVE(P);
	auto TEMP = CURVE.find_nearest(Vertex2D<float32>(1.2,0.99));

	function_bezier2point<float32> ff(CURVE, point);
	newton_minimization_simple(ff, var_arr, LEFT_BORDER, RIGHT_BORDER);
//	Проверка метода нахождения близжайшей точки
	
	TEMP = CURVE.find_nearest(point);
	var_arr.clear();
	x1 = 0.5;
	x2 = 0.5;
	var_arr.push_back(x1);
	var_arr.push_back(x2);
	objective_function<float32> function(CURVE,point);
	newton_minimization<float32>(function, var_arr);
//	*******************************************

    return 0;
	
}

