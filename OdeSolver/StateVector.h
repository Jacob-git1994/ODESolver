#pragma once

#include "OdeSolverParams.h"

#include <valarray>

using std::valarray;
using vec = valarray<double>;

class StateVector
{
private:

	//Our state vector for the current time
	vec currentState;

	//Current params associated with the current state
	OdeSolverParams currentParams;

public:

	//Constructor
	inline StateVector(const valarray<double>&, const OdeSolverParams&);

	//Using default Copy Constructor
	inline StateVector(const StateVector&) = default;

	//Using default Move Constructor
	inline StateVector(StateVector&&) = default;

	//Default Assign operator
	inline StateVector& operator=(const StateVector&) = default;

	//Default Move Assign Operator
	inline StateVector& operator=(StateVector&&) = default;

	//Default Delete operator
	inline ~StateVector() = default;

	//Get the state
	inline const vec& getState() const { return currentState; };
	
	//Get the parameters
	inline const OdeSolverParams& getParams() const { return currentParams; };

};

StateVector::StateVector(const valarray<double>& currentStateIn, const OdeSolverParams& currentParamsIn) :
	currentState(currentStateIn),
	currentParams(currentParamsIn)
{
	//Nothing else to do here
}

