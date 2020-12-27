#pragma once
#include <exception>
#include <iostream>
#include <stdexcept>
#include <memory>
#include <map>
#include <valarray>

#include "SolverIF.h"

//Using these to simplify typing
using std::cerr;
using std::map;
using std::unique_ptr;
using std::exception;
using std::valarray;
using vec = valarray<double>;
using methodPtr = unique_ptr<SolverIF>;
using methodMap = map<unsigned int, methodPtr>;
using std::invalid_argument;

class MethodWrapperIF
{
private:

	//Our method map
	methodMap methods;

	//Build up the solvers
	virtual void buildSolvers() = 0;

protected:

	//Get a handle to the method map
	inline methodMap& getMethodMap() { return methods; };

public:

	//Initalize our method
	virtual void initalize() = 0;

	//Get the solver we want to use
	inline methodPtr& find(SolverIF::SOLVER_TYPES);

	//Update all the methods vectors for new vector size
	inline void updateForVectorSize(const vec&);

};

//Define the find method pointer
methodPtr& MethodWrapperIF::find(SolverIF::SOLVER_TYPES solver)
{
	//Convert solver types to int for the map
	unsigned int solverInt = static_cast<unsigned int>(solver);

	//Get an iterator to method
	methodMap::iterator foundMethod = methods.find(solverInt);

	//Check if method was found
	if (foundMethod == methods.end())
	{
		throw invalid_argument("Invalid Method");
	}
	else
	{
		return foundMethod->second;
	}
}

void MethodWrapperIF::updateForVectorSize(const vec& state)
{
	//Loop through the methods
	for (methodMap::iterator currentMethodItr = methods.begin(); currentMethodItr != methods.end(); ++currentMethodItr)
	{
		//Get the method we are trying to allocate
		methodPtr& currentMethod = currentMethodItr->second;

		try
		{
			currentMethod.get()->initalize(state);
		}
		catch (exception& e)
		{
			cerr << e.what();
		}
	}
}