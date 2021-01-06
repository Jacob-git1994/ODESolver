#pragma once
#include <exception>
#include <iostream>
#include <stdexcept>
#include <memory>
#include <map>
#include <valarray>

#include "Richardson.h"
#include "SolverIF.h"

//Using these to simplify typing
using std::cerr;
using std::map;
using std::unique_ptr;
using std::exception;
using std::valarray;
using vec = valarray<double>;
using methodPtr = unique_ptr<SolverIF>;
using methodMap = map<unsigned int, methodPtr>;
using tableMap = map<unsigned int, Richardson>;
using std::invalid_argument;

class MethodWrapperIF
{
private:

	//Build up the solvers
	virtual void buildSolvers() = 0;

};
