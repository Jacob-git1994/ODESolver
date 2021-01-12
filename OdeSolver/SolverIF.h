#pragma once

#include "OdeFunIF.h"

#include <valarray>

//Convience for writing out methods
using std::valarray;
using vec = valarray<double>;
using crvec = const vec&;
using rvec = vec&;

class SolverIF
{
protected:

	//All instances will have this state
	vec currentState;

	//Update the current state
	inline void updateCurrentState(crvec initalCondition) { currentState = initalCondition; };

	//Method to update the next time step
	inline void updateTimeStep(const double& dt, double& currentTime) { currentTime += dt; };

private:

	//Method for all members to initalze their solving helper vectors
	virtual void initalizeSolverVectors() = 0;

public:

	//Enumerations for the solver types
	enum class SOLVER_TYPES
	{
		EULER				= 10,
		RUNGE_KUTTA_TWO		= 20,
		RUNGE_KUTTA_FOUR	= 30,
		IMPLICT_EULER		= 40,
		CRANK_NICOLSON		= 50
	};

	//Initalize the method's size
	virtual void initalize(crvec) = 0;

	//Get the next time step for rvec
	virtual rvec update(crvec,
						rvec,
						const double&,
						const double&,
						const double&,
						const OdeFunIF*) = 0;

	//Get the power of the error
	virtual const double getErrorOrder() const = 0;

	//Get the current state
	inline crvec getCurrentState() const { return currentState; };
};

