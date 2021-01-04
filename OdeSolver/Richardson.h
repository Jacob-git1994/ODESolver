#pragma once

#include <exception>
#include <iostream>
#include <valarray>

//Some renaming for convience
using std::valarray;
using vec = valarray<double>;
using vecValArray = valarray<vec>;
using mat = valarray<vecValArray>;
using rvec = vec&;
using crvec = const rvec;
using rmat = mat&;
using crmat = const rmat;
using std::exception;
using std::cerr;

class Richardson
{
private:

	//Our Result Matrix
	mat result;

	//The current error vector calculation
	vec currentErrorVector;

	//Our reduction factor of the step size
	double reductionFactor;

	//Our current step size
	double stepSize;

	//Our table size
	unsigned int N;

public:

	//Using default constrcutor
	Richardson() = default;

	//Copy Constructor default
	Richardson(const Richardson&) = default;

	//Use compiler destructor
	~Richardson() = default;

	//Build up our class
	void BuildTables(const size_t);

	//Initalize our steps
	void initalizeSteps(const double&, const double&);

	//Append the result
	void operator()(const size_t, const size_t, crvec);
};

