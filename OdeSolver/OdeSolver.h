#pragma once

#include <chrono>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <vector>

#include "MethodWrapper.h"
#include "MethodWrapperBase.h"
#include "MethodWrapperIF.h"
#include "OdeSolverParams.h"
#include "OdeFunIF.h"
#include "StateVector.h"
#include "SolverIF.h"
#include "Richardson.h"

using std::map;
using std::vector;
using std::unique_ptr;
using std::exception;
using std::runtime_error;
using std::bad_alloc;
using std::cerr;
using std::pow;
using paramMap = map<unsigned int, OdeSolverParams>;
using resultNode = vector<StateVector>;
using results = map<unsigned int, resultNode>;


class OdeSolver
{
private:

	//Our methods
	unique_ptr<MethodWrapperIF> methods;

	//Ode Params for all methods
	OdeSolverParams generalParams;

	//Map of params to mathods
	paramMap params;

	//Map of results in the form of nodes
	results resultMap;

	//Are all the methods explict
	bool isAllExplict() const;

	//Are all the methods implict
	bool isAllImplict() const;

	//Build the richardson tables
	void runMethod(OdeFunIF*,SolverIF*, Richardson&, crvec, rvec, const OdeSolverParams&, const double, const double);

	//Find the best set of parameters for each method
	void findOptimalParams(crvec, OdeFunIF*, const double, const double);

public:

	//Constructor
	OdeSolver(const OdeSolverParams&);

	//Delete the copy constructor
	OdeSolver(const OdeSolver&) = delete;

	//Delete the move operator
	OdeSolver(OdeSolver&&) = delete;

	//Delete the assigment operator
	const OdeSolver& operator=(const OdeSolver&) = delete;

	//Destructor using default
	~OdeSolver() = default;

	void run(OdeFunIF*, crvec,const double, const double);
};

