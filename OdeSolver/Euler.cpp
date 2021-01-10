#include "Euler.h"

//Initalize the Solving Vectors
void Euler::initalizeSolverVectors()
{
	//Get the ref to current method
	Euler& currentMethod = *this;

	//Update the current solver vector size
	currentMethod.k1.resize(currentMethod.getCurrentState().size());
}

//Initalize the vector
void Euler::initalize(crvec initalCondition)
{
	//Get the ref to current method
	Euler& currentMethod = *this;

	//Update the current state
	currentMethod.updateCurrentState(initalCondition);

	//Initalize the Size of the solver vectors
	currentMethod.initalizeSolverVectors();
}

//Update to the next time step
rvec Euler::update(crvec			previousState,
				   rvec				newState,
				   const double&	dt,
				   const double&	tBegin,
				   const double&	tEnd,
	               const OdeFunIF*	functionVector)
{
	//Update the currentState
	currentState = previousState;

	//Save the current time
	double currentTime = tBegin;

	//Iterate through time
	while (currentTime <= tEnd)
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