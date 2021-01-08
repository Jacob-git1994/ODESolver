#pragma once

#include "MethodWrapperIF.h"

//This class will set up the other types of method wrappers
class MethodWrapperBase : public MethodWrapperIF
{
private:

	//Our method map
	methodMap methods;

	//Our table map
	tableMap tables;

	//Build up the solvers
	virtual void buildSolvers() override;

	//Build the tables
	void buildTables();

public:

	//Initalize our method
	void initalize();

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

};
