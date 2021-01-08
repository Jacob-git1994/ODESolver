#include "MethodWrapperBase.h"

//Define the find method pointer
methodPtr& MethodWrapperBase::findMethod(SolverIF::SOLVER_TYPES solver)
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

//Define the find table for the method
Richardson& MethodWrapperBase::findTable(SolverIF::SOLVER_TYPES solver)
{
	//Convert solver types to int for the map
	unsigned int solverInt = static_cast<unsigned int>(solver);

	//Get an iterator to method
	tableMap::iterator foundTable = tables.find(solverInt);

	//Check if method was found
	if (foundTable == tables.end())
	{
		throw invalid_argument("Invalid Method");
	}
	else
	{
		return foundTable->second;
	}
}

void MethodWrapperBase::updateForVectorSize(const vec& state)
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

//Get the map to the methods
methodMap& MethodWrapperBase::getMethodMap()
{
	return methods;
}

//Get the map to the tables
tableMap& MethodWrapperBase::getTableMap()
{
	return tables;
}

//Get the const map to the methods
const methodMap& MethodWrapperBase::getMethodMap() const
{
	return methods;
}

//Get the const map to the tables
const tableMap& MethodWrapperBase::getTableMap() const
{
	return tables;
}

void MethodWrapperBase::buildSolvers()
{
	throw std::runtime_error("Base Method Wrapper does not support building solvers");
}

void MethodWrapperBase::buildTables()
{
	//Go through what methods we are using and add a table class to each one
	for (methodMap::iterator methodIter = methods.begin(); methodIter != methods.end(); ++methodIter)
	{
		tables.emplace(methodIter->first, std::move(Richardson()));
	}
}

void MethodWrapperBase::initalize()
{
	//Build the solvers
	buildSolvers();

	//Build the tables
	buildTables();
}

void MethodWrapperBase::updateForRichardsonTables(const size_t tableSize, const double reductionFactor, const double baseStepSize)
{
	//Go through list of tables and build up the tables
	for (tableMap::iterator richItr = tables.begin(); richItr != tables.end(); ++richItr)
	{
		//Current table
		Richardson& currentTable = richItr->second;

		//Current vector size
		size_t currentVectorSize = findMethod(static_cast<SolverIF::SOLVER_TYPES>(richItr->first)).get()->
			getCurrentState().size();

		//Initalize the current tables parameters
		currentTable.initalizeSteps(reductionFactor, baseStepSize);

		//Initalize the current tables element size
		currentTable.BuildTables(tableSize, currentVectorSize);
	}
}

void MethodWrapperBase::updateAll(const vec& state, const size_t tableSize, const double reductionFactor, const double baseStepSize)
{
	//Update the vector sizes
	updateForVectorSize(state);

	//Update each richardson table
	updateForRichardsonTables(tableSize, reductionFactor, baseStepSize);
}