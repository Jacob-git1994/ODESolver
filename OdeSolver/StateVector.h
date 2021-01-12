#pragma once

#include <valarray>

using std::valarray;
using vec = valarray<double>;

class StateVector
{
private:

	//Our states current time
	double currentTime;

	//Our current error for the result
	double currentError;

	//Our state vector for the current time
	vec currentState;

public:

	//Constructor
	inline StateVector(const double, const double, const vec);

	//Copy Constructor
	inline StateVector(const StateVector&);

	//Assign operator
	inline const StateVector& operator=(const StateVector&);

	//Delete operator using default args
	inline ~StateVector() = default;

	//Get the state
	inline const vec& getState() const;

	//Get the time
	inline const double& getTime() const;

	//Get the error
	inline const double& getError() const;

};

StateVector::StateVector(const double time, const double error, const vec state) :
	currentTime(time),
	currentState(state),
	currentError(error)
{
	//Nothing else to do here
}

StateVector::StateVector(const StateVector& state) :
	currentTime(state.currentTime),
	currentState(state.currentState),
	currentError(state.currentError)
{
	//Nothing else to do here
}

const StateVector& StateVector::operator=(const StateVector& state)
{
	currentTime = state.currentTime;
	currentState = state.currentState;
	currentError = state.currentError;

	return *this;
}

const vec& StateVector::getState() const
{
	return currentState;
}

const double& StateVector::getTime() const
{
	return currentTime;
}

const double& StateVector::getError() const
{
	return currentError;
}