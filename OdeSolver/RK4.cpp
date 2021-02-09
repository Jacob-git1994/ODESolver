#include "RK4.h"

void RK4::initalize(crvec initalCondition)
{
	//Call Base Classes Initalization
	Euler::initalize(initalCondition);
}

void RK4::initalizeSolverVectors()
{
	//Get a handle to this
	RK4& currentMethod = *this;

	//Initalize our solvers
	currentMethod.k1.resize(currentMethod.currentState.size());
	currentMethod.k2.resize(currentMethod.currentState.size());
	currentMethod.k3.resize(currentMethod.currentState.size());
	currentMethod.k4.resize(currentMethod.currentState.size());
}

rvec RK4::update(crvec previousState, rvec newState, const double& dt, const double& beginTime, const int& numOfSteps, const OdeFunIF* problem)
{
	//Update the current state
	currentState = previousState;

	//Save the current time
	double currentTime = beginTime;

	//Iterate through time
	for (int i = 0; i < numOfSteps; ++i)
	{
		//Update all our solvers
		k1 = problem->operator()(k1, currentState, currentTime);
		k2 = problem->operator()(k2, currentState + dt * (k1 / 2.0), currentTime + dt / 2.0);
		k3 = problem->operator()(k3, currentState + dt * (k2 / 2.0), currentTime + dt / 2.0);
		k4 = problem->operator()(k4, currentState + dt * k3, currentTime + dt);

		//Update the current state with the weighted average
		currentState += (dt / 6.) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);

		//Update the time step
		updateTimeStep(dt, currentTime);
	}

	//Save the currentState
	newState = currentState;

	//Return the new state
	return newState;
}

/// <summary>
/// Will run the implict solver. Will fail as RK4 is not implict.
/// </summary>
/// <param name="previousState"></param>
/// <param name="newState"></param>
/// <param name="dt"></param>
/// <param name="beginTime"></param>
/// <param name="numOfSteps"></param>
/// <param name="problem"></param>
/// <param name="implictDt"></param>
/// <param name="implictError"></param>
/// <returns></returns>
rvec RK4::update(crvec previousState, rvec newState, const double& dt, const double& beginTime, const int& numOfSteps, const OdeFunIF* problem, const double& implictDt, const double& implictError)
{
	throw logic_error("Implict Method Not Implimented in Explict Scheme");
}

const double RK4::getErrorOrder() const
{
	return 4.0;
}