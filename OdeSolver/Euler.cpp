#include "Euler.h"

/// <summary>
/// Initalize the vectors (k1,k2...) to be used to solve this system
/// </summary>
void Euler::initalizeSolverVectors()
{
	//Get the ref to current method
	Euler& currentMethod = *this;

	//Update the current solver vector size
	currentMethod.k1.resize(currentMethod.getCurrentState().size());
}

/// <summary>
/// Update the current state with the inital condition so we know the size and update the solving helper vectors (funcs evaled at different "future" times)
/// </summary>
/// <param name="initalCondition"></param>
void Euler::initalize(crvec initalCondition)
{
	//Get the ref to current method
	Euler& currentMethod = *this;

	//Update the current state
	currentMethod.updateCurrentState(initalCondition);

	//Initalize the Size of the solver vectors
	currentMethod.initalizeSolverVectors();
}

/// <summary>
/// Find the state vector at the next time step defined by the function derrivative vector
/// </summary>
/// <param name="previousState"></param>
/// <param name="newState"></param>
/// <param name="dt"></param>
/// <param name="tBegin"></param>
/// <param name="numOfSteps"></param>
/// <param name="functionVector"></param>
/// <returns></returns>
rvec Euler::update(crvec			previousState,
				   rvec				newState,
				   const double&	dt,
				   const double&	tBegin,
				   const int&		numOfSteps,
	               const OdeFunIF*	functionVector)
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

		//Update to the next time step and save the current state
		currentState += dt * k1;

		//Update the time step to the next time
		updateTimeStep(dt, currentTime);
	}

	//Save off the final current state to the new state
	newState = currentState;

	//Return the new state
	return newState;
}

/// <summary>
/// This runs the implict calculations. We will throw here as euler is not implict
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
rvec Euler::update(crvec previousState, rvec newState, const double& dt, const double& beginTime, const int& numOfSteps, const OdeFunIF* problem, const double& implictDt, const double& implictError)
{
	throw logic_error("Implict Method Not Implimented in Explict Scheme");
}

/// <summary>
/// Get the leading error order
/// </summary>
/// <returns></returns>
const double Euler::getErrorOrder() const
{
	return 2.;
}