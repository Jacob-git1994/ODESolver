#pragma once

#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <thread>
#include <vector>

#include "MethodWrapperBase.h"
#include "OdeSolverParams.h"
#include "OdeFunIF.h"
#include "StateVector.h"
#include "SolverIF.h"
#include "Richardson.h"

using std::map;
using std::vector;
using std::shared_ptr;
using std::exception;
using std::runtime_error;
using std::bad_alloc;
using std::cerr;
using std::pow;
using std::thread;
using std::setprecision;
using paramMap = map<unsigned int, OdeSolverParams>;
using resultNode = vector<StateVector>;
using results = map<unsigned int, resultNode>;
using threadVector = vector<thread>;


class OdeSolver
{
private:

	//Our methods
	MethodWrapperBase methods;

	//Ode Params for all methods
	OdeSolverParams generalParams;

	//Map of params to mathods
	paramMap params;

	//Map of results in the form of nodes
	results resultMap;

	//vector to store threads
	threadVector methodThreads;

	//Set up all our methods and tables
	void setup();

	//Build the richardson tables
	void runMethod(const OdeFunIF*, unique_ptr<SolverIF>&, Richardson&, crvec, rvec, const OdeSolverParams&, const double, const double);

	//Once the best parameters are founds we can start solving for the next time step
	vec buildSolution(unique_ptr<SolverIF>&, const unsigned int, crvec, const OdeFunIF*, const double, const double);

	//Update the tables
	void updateMethod(unique_ptr<SolverIF>&, const OdeSolverParams&, Richardson&, crvec, rvec, const double, const double, const double, const OdeFunIF*, const int);

	//Update dt
	const bool updateDt(OdeSolverParams&, const bool, const double, const double);

public:

	//Delete the default constructor
	OdeSolver() = delete;

	//Constructor
	OdeSolver(const OdeSolverParams&);

	//Delete the copy constructor
	OdeSolver(const OdeSolver&) = delete;

	//Delete the assignment operator
	OdeSolver& operator=(const OdeSolver&) = delete;

	//Use the default move constructor
	OdeSolver(OdeSolver&&) = default;

	//Use the default move assigment constructor
	OdeSolver& operator=(OdeSolver&&) = delete;

	//Destructor using default
	~OdeSolver() = default;

	//Run our method
	void run(OdeFunIF*, crvec,const double, const double, const unsigned int);

	//Clear out our data for another run
	void refreshParams(const OdeSolverParams&);
};

