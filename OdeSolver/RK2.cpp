#include "RK2.h"

/// <summary>
/// Initalize our current state vector - call previous method Euler as functionality does not change
/// </summary>
/// <param name="initalConditionsIn"></param>
void RK2::initalize(const valarray<double>& initalConditionsIn)
{
	Euler::initalize(initalConditionsIn);
}

void RK2::initalizeSolverVectors()
{
	//Get ref to our current method
	RK2& currentMethod = *this;

	//Update the current solver vector size
	currentMethod.k1.resize(currentMethod.getCurrentState().size());
	currentMethod.k2.resize(currentMethod.getCurrentState().size());
}

valarray<double>& RK2::update(
	const valarray<double>& previousState,
	valarray<double>& newState,
	const double& dt,
	const double& tBegin,
	const int& numOfSteps,
	const OdeFunIF* functionVector)
{
	//Update the currentState
	currentState = previousState;

	//Save the current time
	double currentTime = tBegin;

	//Iterate through time
	for (int i = 0; i < numOfSteps; ++i)
	{
		//Get the function vector at the current time
		k1 = functionVector->operator()(k1, currentState, currentTime);
		k2 = functionVector->operator()(k2, currentState + 0.5 * dt * k1, currentTime + 0.5 * dt);

		//Update to the next time step and save the current state
		currentState += dt * k2;

		//Update the time step to the next time
		updateTimeStep(dt, currentTime);
	}

	//Save off the final current state to the new state
	newState = currentState;

	//Return the new state
	return newState;
}

const double RK2::getErrorOrder() const
{
	return 3.0;
}