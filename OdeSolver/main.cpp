#include "Euler.h"
#include "OdeFunIF.h"
#include "SolverIF.h"
#include "MethodWrapper.h"
#include "Richardson.h"

#include <valarray>
#include <iostream>

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
	state[0] = 1. * currentState[0];

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

	MethodWrapper methods;
	methods.initalize();
	methods.updateForVectorSize(ic);

	Richardson testCleaning;
	testCleaning.initalizeSteps(2., .01);
	testCleaning.BuildTables(10,ic.size());
	

	for (int i = 0; i < 10; ++i)
	{
		methods.updateForVectorSize(ic);
		methods.getSolver()->update(sol, .01/pow(2.,i), initalTime, 1., testProblem);
		//std::cout << sol[0] << "\n";
		testCleaning(i, 0, sol);
	}

	std::valarray<double> good;
	cout << testCleaning.error(good) << "\n";
	cout << good[0] << "\n";

	delete testProblem;
}