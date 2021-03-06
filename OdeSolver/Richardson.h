#pragma once

#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <valarray>

#include "SolverIF.h"

//Some renaming for convience
using vecValArray = valarray<vec>;
using mat = valarray<vecValArray>;
using rmat = mat&;
using crmat = const rmat;
using std::exception;
using std::cerr;
using std::sqrt;
using std::pow;
using std::fabs;
using std::log;
using std::invalid_argument;
using std::for_each;

class Richardson
{
private:

	//Calculate the vector norm
	double normedError() const;

	//Our Result Matrix
	mat result;

	//The current error vector calculation
	double currentNormError = 0.0;

	//Our reduction factor of the step size
	double reductionFactor = 0.0;

	//Our current step size
	double stepSize = 0;

	//Our table size
	unsigned int N = 0;

	//Flag to check if tables are built
	bool isBuilt = false;

public:

	//Using default constrcutor (need to add something here to initalize parameters)
	Richardson() = default;

	//Copy Constructor default
	Richardson(const Richardson&) = default;

	//Use compiler destructor
	~Richardson() = default;

	//Build up our class
	void BuildTables(const size_t, const size_t);

	//Initalize our steps
	void initalizeSteps(const double&, const double&);

	//Append result moving the result
	void append(const size_t, const size_t, valarray<double>&&);

	//Append the result copying the result
	void append(const size_t, const size_t, const valarray<double>&);

	//Get the error, updated vector, and estimate of the orders constant
	const double error(rvec, double&);

	//Get the table size
	const size_t getTableSize() const;

	//Get the reduction factor
	const double getReductionFactor() const;
};

