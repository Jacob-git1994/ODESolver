#pragma once

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
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
using std::exception;
using std::runtime_error;
using std::bad_alloc;
using std::unique_ptr;
using std::cerr;
using std::pow;
using std::thread;
using std::setprecision;
using paramMap = map<unsigned int, OdeSolverParams>;
using resultNode = vector<StateVector>;
using results = map<unsigned int, resultNode>;
using threadVector = vector<thread>;

// Class to hold all the methods, results, and parameters.
// This class will also build up all the methods and distribute the methods to find a new solution on each parallel thread.
// Post-Run, this class will support retreaving the data generated by the whole array of states or an approximation to a paticular value.
class OdeSolver
{
private:

	// The methods class will hold all the methods needed and their richardson tables to be run. 
	MethodWrapperBase methods;

	// This is the base set of parameters the user has selected to be run.
	// These will be distributed to each method by Id in params
	OdeSolverParams generalParams;

	// This is the map of all the parameters to each method. Each run we will save off the parameters generated by the method
	paramMap params;

	// This is the map of the entire solution to the current method. 
	// The state vector will hold the state generated at the time and parameters used to solve for that time. 
	map<unsigned int, vector<StateVector>> resultMap;

	// This is a vector of threads that we will use to solve the problem in parallell for each method. 
	// These will be updated in the core running section.
	threadVector methodThreads;

	// This starts up saving all the parameters and seeing which methods the user wants.
	// It will call on methodbasewrapper to build each method of what is allowed and build each method with a corresponding richardson table.
	// It will also generate the result map and parameter map for each allowable method.
	void setup();

	// This will run the paticular method referenced in input arguments.
	void runMethod(const OdeFunIF*, unique_ptr<SolverIF>&, Richardson&, crvec, rvec, const OdeSolverParams&, const double, const double);

	// This will build the solution from the current time step to the next "best" time step.
	vec buildSolution(unique_ptr<SolverIF>&, const unsigned int, Richardson&, OdeSolverParams&, crvec, const OdeFunIF*, const double, const double);

	// This updates the method to the next time step. 
	// This is used in each thread. 
	void updateNextTimeStep(
		const unsigned int,
		unique_ptr<SolverIF>&,
		OdeSolverParams&,
		Richardson&,
		const double,
		const double,
		const valarray<double>&,
		const OdeFunIF*,
		vector<StateVector>&);

	// This method calls the current method's method to update the current state to the next time step. 
	void updateMethod(unique_ptr<SolverIF>&, const OdeSolverParams&, Richardson&, crvec, rvec, const double, const double, const double, const OdeFunIF*, const int);

	// Check the error and determine if an upgrade or downgrade is required to satify the current estimated error. 
	// If we fail and we are not on the last iteration, we will find the new dt and run the iteration scheme again
	const bool updateDt(OdeSolverParams&, const bool, const double, const double);

public:

	//Delete the default constructor
	OdeSolver() = default;

	//Constructor
	OdeSolver(const OdeSolverParams&);

	//Delete the copy constructor
	OdeSolver(const OdeSolver&) = delete;

	//Delete the assignment operator
	OdeSolver& operator=(const OdeSolver&) = delete;

	//Use the default move constructor
	OdeSolver(OdeSolver&&) noexcept = default;

	//Use the default move assigment constructor
	OdeSolver& operator=(OdeSolver&&) noexcept = default;

	//Destructor using default
	~OdeSolver() = default;

	//Run our method
	void run(const OdeFunIF*, crvec, const double, const double);

	//Clear out our data for another run
	void refreshParams(const OdeSolverParams&);

	//Get the results for a given type
	const vector<StateVector>& getResults(SolverIF::SOLVER_TYPES) const;

	//Get the results with an unsigned int if the enums are known
	const vector<StateVector>& getResults(const unsigned int) const;

	//Get the results of the best method
	const vector<StateVector>& getResults() const;

	//Get the results with a paticular method and time
	StateVector getStateAndTime(SolverIF::SOLVER_TYPES, const double) const;

	//Get the results with a paticular method known enum value and time
	StateVector getStateAndTime(const unsigned int, const double) const;

	//Find the best result and return the state vector interplation of that result
	StateVector getStateAndTime(const double) const;
};

