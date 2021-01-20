#include "MethodRK4.h"

void MethodRK4::buildSolvers()
{
	try
	{
		//Add Euler method
		getMethodMap().emplace(static_cast<unsigned int>(SolverIF::SOLVER_TYPES::RUNGE_KUTTA_FOUR),
			std::move(unique_ptr<SolverIF>(new RK4)));
	}
	catch (exception& e)
	{
		cerr << e.what();
	}
}