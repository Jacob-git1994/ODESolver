#pragma once
#include <valarray>

//Convience for writing out methods
using std::valarray;
using vec = valarray<double>;
using crvec = const vec&;
using rvec = vec&;

class OdeFunIF
{
public:

	//Just need to define this operator to pass to a new ode class
	virtual rvec operator()(rvec,
							crvec,
							const double&) const = 0;
};

