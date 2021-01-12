#pragma once

#include <cmath>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <valarray>

#include "SolverIF.h"

//Some renaming for convience
//using std::valarray;
//using vec = valarray<double>;
using vecValArray = valarray<vec>;
using mat = valarray<vecValArray>;
//using rvec = vec&;
//using crvec = const rvec;
using rmat = mat&;
using crmat = const rmat;
using std::exception;
using std::cerr;
using std::sqrt;
using std::pow;
using std::fabs;
using std::log;
using std::invalid_argument;

class Richardson
{
private:

	//Our Result Matrix
	mat result;

	//The current error vector calculation
	double currentNormError;

	//Our reduction factor of the step size
	double reductionFactor;

	//Our current step size
	double stepSize;

	//Our table size
	unsigned int N;

	//Calculate the vector norm
	double normedError() const;

	//Flag to check if tables are built
	bool isBuilt;

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

	//Append the result
	void operator()(const size_t, const size_t, crvec);

	//Get the error, updated vector, and estimate of the orders constant
	double error(rvec, double&, const double);

	//Get the error ignoring the updated vector
	double error(double&, const double);

	//Just get the error founds ignorning the other elements
	double error();

	//Get the table size
	const size_t getTableSize() const;

	//Get the reduction factor
	const double getReductionFactor() const;
};

