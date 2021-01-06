#include "MethodWrapper.h"

//Build up our solvers (Start with just building one for right now
void MethodWrapper::buildSolvers()
{
	//Check to make sure we allocated properly
	try
	{
		//Add Euler method
		getMethodMap().emplace(static_cast<unsigned int>(SolverIF::SOLVER_TYPES::EULER),
							   std::move(unique_ptr<SolverIF>(new Euler)));
	}
	catch (exception& e)
	{
		cerr << e.what();
	}
	//If this fails move onto the next method
}