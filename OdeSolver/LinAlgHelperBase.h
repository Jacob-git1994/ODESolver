#pragma once

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <valarray>

#include "OdeFunIF.h"

using std::valarray;

class LinAlgHelperBase
{
protected:

	//Matrix A
	valarray<valarray<double>> A;

	//Our left guess
	valarray<double> guessLeft;

	//Our right guess
	valarray<double> guessRight;

	//Our function vector
	valarray<double> funcVec;

	//Error tolerance
	double errorTol;

	//Max number of iterations
	unsigned int maxIter;

	//For the partial derivatives
	double dt;

	//Generate the Jacobian
	virtual void getJacobian(const OdeFunIF*, const double&, const double&) = 0;

	//Generate the function derivative vector
	void getFuncDer(const OdeFunIF*, const double&);

	//Solver the system
	const valarray<double> solveSystem();

public:

	//Constructor del
	LinAlgHelperBase() = delete;

	//This will be generated when constructing each method for each implict method
	LinAlgHelperBase(const double&, const double&, const unsigned int&);

	//Default copy constructor
	LinAlgHelperBase(const LinAlgHelperBase&) = default;

	//Default assin constructor
	LinAlgHelperBase& operator=(const LinAlgHelperBase&) = default;

	//Default destructor
	virtual ~LinAlgHelperBase() = default;

	//Solve the problem
	const valarray<double>& solve(const double&, const double&, const valarray<double>&, const OdeFunIF*);
	
};

