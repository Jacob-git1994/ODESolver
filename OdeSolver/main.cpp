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
	//Generate velocity componets
	state[0] = currentState[3];
	state[1] = currentState[4];
	state[2] = currentState[5];

	//Generate acceleration componets
	state[3] = 0;
	state[4] = 0;
	state[5] = -(5.972e+27 * 10. * 6.673e-11) / (currentState[2] * currentState[2] + 6.3781e6 * 6.3781e6);

	//Reset the state if we go negetive
	for (int i = 0; i < 3; ++i)
	{
		//Hit the earth
		if (state[i] < 0)
		{
			state[i] = 0;
			state[i + 3] = 0.0;
		}
	}

	return state;
}

int main()
{
	//Initalize our test problem
	OdeFunIF* testProblem = new Test;

	//Initalize our inital condion
	std::valarray<double> ic = { 0,0,0,0,0,10000};
	std::valarray<double> sol(6);
	double initalTime = 0.;

	/*
	OdeSolverParams params({ true,false,false,false,false },
		{ .001,.01 },
		{ 10,60 },
		{ 2,8 },
		{ false,false,false });

	*/
	OdeSolverParams params;

	params.upperError = 1e-5;
	params.lowerError = 1e-6;
	params.redutionFactor = 2.;
	params.dt = .01;
	params.minDt = .1;
	params.maxDt = 10.;
	params.minTableSize = 4;
	params.maxTableSize = 8;
	params.useEuler = false;
	params.useRK4 = true;
	params.useRK2 = false;
	params.smallestAllowableDt = 1e-4;

	OdeSolver solver(params);
	OdeSolver solv2;

	solv2 = std::move(params);

	solver.refreshParams(params);
		
	solver.run(testProblem, ic, 0.0, 1000);

	const unsigned int maxNodes = 1000;
	for (int i = 0; i <= maxNodes; ++i)
	{
		const double step = static_cast<double>(i);// / static_cast<double>(maxNodes);
		std::cout << solver.getStateAndTime(step).getParams().currentTime << "\t" << solver.getStateAndTime(step).getState()[2] << "\t" 
			<< solver.getStateAndTime(step).getParams().totalError << "\t" <<  solver.getResults().size()  << "\n";
	}

	//std::cout << solver.getStateAndTime(1).getState()[0] << "\t" << solver.getStateAndTime(1).getParams().totalError << "\n";

	delete testProblem;
}