#pragma once

#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <valarray>

#include "Euler.h"
#include "Richardson.h"
#include "RK4.h"
#include "OdeSolverParams.h"

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

//This class will set up the other types of method wrappers
class MethodWrapperBase
{
private:

	//Build up the solvers
	void buildSolvers(const OdeSolverParams&);

	//Build the tables
	void buildTables();

	//Our method map
	methodMap methods;

	//Our table map
	tableMap tables;

public:

	//Using default constructor
	MethodWrapperBase() = default;

	//Delete the copy constructor
	MethodWrapperBase(const MethodWrapperBase&) = delete;

	//Delete the assignment operator
	MethodWrapperBase& operator=(const MethodWrapperBase&) = delete;

	//Using default move operator
	MethodWrapperBase(MethodWrapperBase&&) = default;

	//Using default assign move operator
	MethodWrapperBase& operator=(MethodWrapperBase&&) = default;

	//Default destructor
	~MethodWrapperBase() = default;

	//Initalize our method
	void initalize(const OdeSolverParams&);

	//Get a referance to the method map
	methodMap& getMethodMap();

	//Get a referance to the table map
	tableMap& getTableMap();

	//Get a const referance to the method map
	const methodMap& getMethodMap() const;

	//Get a const referance to the table map
	const tableMap& getTableMap() const;

	//Get the solver we want to use
	methodPtr& findMethod(SolverIF::SOLVER_TYPES);

	//Get pointer to tables
	Richardson& findTable(SolverIF::SOLVER_TYPES);

	//Update all the methods vectors for new vector size
	void updateForVectorSize(const vec&);

	//Update all tables for the Richardson Table Sizes
	void updateForRichardsonTables(const size_t, const double, const double);

	//Update all the methods vectors for new vector size
	void updateAll(const vec&, const size_t, const double, const double);

	//Clear out all the methods we used
	void clearMethods();
};
