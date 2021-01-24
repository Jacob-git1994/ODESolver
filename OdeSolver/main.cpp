#include "Euler.h"
#include "OdeFunIF.h"
#include "SolverIF.h"
#include "Richardson.h"
#include "OdeSolverParams.h"
#include "OdeSolver.h"

#include <valarray>
#include <iostream>
#include <cmath>

//using vec = std::valarray<double>;
using std::cout;

//Create a test solver
class Test : public OdeFunIF
{
public:

	virtual std::valarray<double>& operator()(std::valarray<double>&,
							const std::valarray<double>&,
							const double&) const override;
	
};

//Define our updating method
std::valarray<double>& Test::operator()(std::valarray<double>& state,
					  const std::valarray<double>& currentState,
					  const double& currentTime) const
{
	state[0] = currentState[0];

	return state;
}

int main()
{
	//Initalize our test problem
	OdeFunIF* testProblem = new Test;

	//Initalize our inital condion
	std::valarray<double> ic = {1.};
	std::valarray<double> sol(1);
	double initalTime = 0.;

	/*
	OdeSolverParams params({ true,false,false,false,false },
		{ .001,.01 },
		{ 10,60 },
		{ 2,8 },
		{ false,false,false });

	*/
	OdeSolverParams params;

	params.upperError = 1e-4;
	params.lowerError = 1e-5;
	params.redutionFactor = 2.;
	params.dt = .01;
	params.minDt = .01;
	params.maxDt = 10.;
	params.minTableSize = 8;
	params.maxTableSize = 15;
	params.useEuler = true;
	params.useRK4 = true;
	params.useRK2 = true;

	OdeSolver solver(params);
	OdeSolver solv2;

	solv2 = std::move(params);

	solver.refreshParams(params);
		
	solver.run(testProblem, ic, 0.0, 1,1000);

	std::cout << solver.getStateAndTime(SolverIF::SOLVER_TYPES::RUNGE_KUTTA_FOUR, 1).getState()[0] << "\t" << solver.getStateAndTime(SolverIF::SOLVER_TYPES::RUNGE_KUTTA_FOUR, 1).getParams().totalError << "\n";

	delete testProblem;
}