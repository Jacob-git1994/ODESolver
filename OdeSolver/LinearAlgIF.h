#pragma once

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <valarray>

#include "OdeFunIF.h"

using std::valarray;
using std::for_each;

class LinearAlgIF
{
protected:

	//Storage for the Jacobian
	valarray<valarray<double>> J;

	//Our left result
	valarray<double> leftResult;

	//Our right result
	valarray<double> rightResult;

	//Our left function derivative
	valarray<double> leftDfDt;

	//Our right function derivative
	valarray<double> rightDfDt;

	//Control the iterations
	unsigned int maxIterAllowed;

	//Error Tolerance
	double allowedError;

	//Our step size
	double stepDt;

private:

	//Build our Jacobian
	virtual void buildJacbian(const OdeFunIF*, const valarray<double>&, const double&) = 0;

	//Build our function vector
	void buildFuncDer(const OdeFunIF*, valarray<double>&, const valarray<double>&, const double&);

	//Solve our system of linear equations
	const valarray<double> solveSystem();

public:

	//Default constructor
	LinearAlgIF() = default;

	//Default copy constructor
	LinearAlgIF(const LinearAlgIF&) = default;

	//Default assign constructor
	LinearAlgIF& operator=(const LinearAlgIF&) = default;

	//Initalize our variables (vector size & maxiterallowed)
	void initalize(const size_t&, const double&);
};

