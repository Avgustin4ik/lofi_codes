#pragma once
template <typename F,typename T>

class derivative_r //(f(x + h) - f(x)) / h
{
public:
	derivative_r(const F& f, const T& h) :f(f), h(h) {}
	Vertex2D<T> operator () (const T& x)
	{
		return (f(x + h) - f(x)) / h;
	}
	const Vertex2D<T> operator () (const T& x) const
	{
		return (f(x + h) - f(x)) / h;
	}
private:
	const F& f;
	T h;
};

template <typename F, typename T>

class derivative_l //(f(x) - f(x - h)) / h
{
public:
	derivative_l(const F& f, const T& h) :f(f), h(h) {}
	Vertex2D<T> operator () (const T& x)
	{
		return (f(x) - f(x - h)) / h;
	}
	const Vertex2D<T> operator () (const T& x) const
	{
		return (f(x) - f(x - h)) / h;
	}
private:
	const F& f;
	T h;
};

template <typename F, typename T>

class derivarive_2p //производная по двум точкам (слева и права) (f(x + h) - f(x - h)) / (2 * h)
{
public:
	derivarive_2p(const F& f, const T& h) :f(f), h(h) {}
	Vertex2D<T> operator () (const T& x)
	{
		return (f(x + h) - f(x - h)) / (float(2) * h);
	}
	const Vertex2D<T> operator () (const T& x) const
	{
		return (f(x + h) - f(x - h)) / (float(2) * h);
	}
private:
	const F& f;
	T h;
};
