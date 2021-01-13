#include "Euler.h"
#include "OdeFunIF.h"
#include "SolverIF.h"
#include "MethodWrapper.h"
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
	state[0] = std::exp(currentTime);//std::cos(currentTime);

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

	params.upperError = .0000009;
	params.lowerError = .00000001;
	params.redutionFactor = 4.;
	params.dt = .01;
	params.minDt = .00009;
	params.maxDt = .0001;
	params.minTableSize = 4;
	params.maxTableSize = 6;

	OdeSolver solver(params);
		
	solver.run(testProblem, ic, 0, 1.);

	delete testProblem;
}