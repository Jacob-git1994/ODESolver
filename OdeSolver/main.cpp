#include "Euler.h"
#include "OdeFunIF.h"
#include "SolverIF.h"
#include "MethodWrapper.h"

#include <valarray>
#include <iostream>

using vec = std::valarray<double>;
using std::cout;

//Create a test solver
class Test : public OdeFunIF
{
public:

	virtual vec& operator()(vec&,
							const vec&,
							const double&) const override;
	
};

//Define our updating method
vec& Test::operator()(vec& state,
					  const vec& currentState,
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
	vec ic = {1.};
	vec sol(1);
	double initalTime = 0.;

	MethodWrapper methods;
	methods.initalize();
	methods.updateForVectorSize(ic);
	methods.getSolver()->update(sol, .000001, initalTime, 1., testProblem);

	cout << methods.getSolver()->getCurrentState()[0] << "\n";

	delete testProblem;
}