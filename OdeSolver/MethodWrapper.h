#pragma once
#include <exception>
#include <iostream>
#include <map>
#include <memory>
#include <valarray>

#include "Euler.h"
#include "MethodWrapperIF.h"
#include "SolverIF.h"

//Using these to simplify typing
using std::exception;
using std::map;
using std::unique_ptr;
using std::make_unique;
using std::cerr;
using methodPtr = unique_ptr<SolverIF>;
using methodMap = map<unsigned int, methodPtr>;

class MethodWrapper : public MethodWrapperIF
{
private:

	//build our solvers for these types of methods
	virtual void buildSolvers() override;

public:

	//Default Constructor
	MethodWrapper() = default;

	//Default Copy consutrctor
	MethodWrapper(const MethodWrapper&) = default;

	//Default deconstructor
	~MethodWrapper() = default;

	//Initalize the solvers
	virtual void initalize() override;

	//Just for a test
	inline SolverIF* getSolver() { return find(SolverIF::SOLVER_TYPES::EULER).get(); };
};

