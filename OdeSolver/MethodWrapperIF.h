#pragma once
#include <exception>
#include <iostream>
#include <stdexcept>
#include <memory>
#include <map>
#include <valarray>

#include "Richardson.h"
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
using tableMap = map<unsigned int, Richardson>;
using std::invalid_argument;

class MethodWrapperIF
{
protected:

	//Build up the solvers
	virtual void buildSolvers() = 0;

public:

	//Initalize our method
	virtual void initalize() = 0;

	//Get a referance to the method map
	virtual methodMap& getMethodMap() = 0;

	//Get a referance to the table map
	virtual tableMap& getTableMap() = 0;

	//Get a const referance to the method map
	virtual const methodMap& getMethodMap() const = 0;

	//Get a const referance to the table map
	virtual const tableMap& getTableMap() const = 0;

	//Get the solver we want to use
	virtual methodPtr& findMethod(SolverIF::SOLVER_TYPES) = 0;

	//Get pointer to tables
	virtual Richardson& findTable(SolverIF::SOLVER_TYPES) = 0;

	//Update all the methods vectors for new vector size
	virtual void updateForVectorSize(const vec&) = 0;

	//Update all tables for the Richardson Table Sizes
	virtual void updateForRichardsonTables(const size_t, const double, const double) = 0;

	//Update all the methods vectors for new vector size
	virtual void updateAll(const vec&, const size_t, const double, const double) = 0;

};
