#include "Euler.h"
#include "OdeFunIF.h"
#include "SolverIF.h"
#include "Richardson.h"
#include "OdeSolverParams.h"
#include "OdeSolver.h"

#include <valarray>
#include <iostream>
#include <cmath>
#include <iomanip>

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
	/*
	double thrust = 0.0;
	if (currentTime <= 40)
	{
		thrust = 100;
	}

	double theta = .5;

	state[0] = currentState[2];
	state[1] = currentState[3];
	state[2] = -thrust * std::sin(theta);
	state[3] = thrust * std::cos(theta);

	return state;
	*/
	state[0] = currentState[0];//currentState[0]*std::cos(currentTime)*currentTime;
	return state;
}

int main()
{
	//Initalize our test problem
	OdeFunIF* testProblem = new Test;

	//Initalize our inital condion
	std::valarray<double> ic = {1};//{ 0,0,0,0,0,10000};
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

	params.upperError = 1e-7;
	params.lowerError = 1e-11;
	params.redutionFactor = 2.;
	params.dt = .001;
	params.minDt = .1;
	params.maxDt = 2.;
	params.minTableSize = 15;
	params.maxTableSize = 16;
	params.useEuler = true;
	params.useRK4 = true;
	params.useRK2 = true;
	params.isFast = false;
	params.smallestAllowableDt = 1e-4;

	auto begin = 0.0;
	auto end = 8.0;
	auto dt = .01;

	OdeSolver solver(params);
	OdeSolver solv2 = std::move(params);

	solver.refreshParams(params);
		
	solver.run(testProblem, ic, begin, end);

	const auto grid = [&]()
	{
		std::vector<double> timeGrid;

		auto& currentTime = begin;
		timeGrid.push_back(currentTime);
		for (; currentTime <= end; currentTime += dt)
		{
			timeGrid.push_back(currentTime);
		}
		timeGrid.push_back(end);

		return timeGrid;
	};

	for (auto&& time : grid())
	{
		const auto& state = solver.getStateAndTime(time);

		std::cout << std::setprecision(14) << state.getParams().currentTime << "\t" << state.getState()[0] << std::endl;
	}

	/*
	for (const auto& sol : solver.getResults())
	{
		std::cout << std::setprecision(18) << std::setw(18) << std::left << sol.getParams().currentTime << "\t" << std::left << sol.getState()[0] <<  "\t" << std::left << sol.getParams().dt << "\t" << std::left << sol.getParams().currentTableSize << "\t" << std::left << sol.getParams().totalError << "\t" << std::left << sol.getParams().c << "\n";

		//std::cout << std::setprecision(18) << std::setw(18) << std::left << sol.getParams().currentTime << "\t" << std::left << sol.getState()[0] << "\t" << sol.getState()[1]  << "\t" << sol.getState()[2] << "\t" << sol.getState()[3] << "\t" << std::left << sol.getParams().dt << "\t" << std::left << sol.getParams().currentTableSize<< "\t" << std::left << sol.getParams().totalError << "\t" << std::left << sol.getParams().c << "\n";
	}
	*/

	std::cout << solver.getStateAndTime(.5).getState()[0] << "\n";

	/*
	const unsigned int maxNodes = 1000;
	for (int i = 0; i <= maxNodes; ++i)
	{
		const double step = static_cast<double>(i);// / static_cast<double>(maxNodes);
		std::cout << solver.getStateAndTime(step).getParams().currentTime << "\t" << solver.getStateAndTime(step).getState()[2] << "\t" 
			<< solver.getStateAndTime(step).getParams().totalError << "\t" <<  solver.getResults().size()  << "\n";
	}
	*/

	//std::cout << solver.getStateAndTime(1).getState()[0] << "\t" << solver.getStateAndTime(1).getParams().totalError << "\n";

	delete testProblem;
}