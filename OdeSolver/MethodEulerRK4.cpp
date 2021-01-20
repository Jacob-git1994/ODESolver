#include "MethodEulerRK4.h"

void MethodEulerRK4::buildSolvers()
{
	//Build Euler Solvers
	MethodEuler::buildSolvers();

	//Build runge kutta 4 solvers
	MethodRK4::buildSolvers();
}