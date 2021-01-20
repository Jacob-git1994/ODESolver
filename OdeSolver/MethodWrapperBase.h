#pragma once

#include "MethodWrapperIF.h"

//This class will set up the other types of method wrappers
class MethodWrapperBase : virtual public MethodWrapperIF
{
protected:

	//Build up the solvers
	virtual void buildSolvers() override;

private:

	//Our method map
	methodMap methods;

	//Our table map
	tableMap tables;

	//Build the tables
	void buildTables();

public:

	//Using default constructor
	MethodWrapperBase() = default;

	//Delete the copy constructor
	MethodWrapperBase(const MethodWrapperBase&) = delete;

	//Default destructor
	virtual ~MethodWrapperBase() = default;

	//Initalize our method
	virtual void initalize() override;

	//Get a referance to the method map
	virtual methodMap& getMethodMap() override;

	//Get a referance to the table map
	virtual tableMap& getTableMap() override;

	//Get a const referance to the method map
	virtual const methodMap& getMethodMap() const override;

	//Get a const referance to the table map
	virtual const tableMap& getTableMap() const override;

	//Get the solver we want to use
	virtual methodPtr& findMethod(SolverIF::SOLVER_TYPES) override;

	//Get pointer to tables
	virtual Richardson& findTable(SolverIF::SOLVER_TYPES) override;

	//Update all the methods vectors for new vector size
	virtual void updateForVectorSize(const vec&) override;

	//Update all tables for the Richardson Table Sizes
	virtual void updateForRichardsonTables(const size_t, const double, const double) override;

	//Update all the methods vectors for new vector size
	virtual void updateAll(const vec&, const size_t, const double, const double) override;

};
