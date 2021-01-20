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

void MethodWrapperBase::buildTables()
{
	//Go through what methods we are using and add a table class to each one
	for (methodMap::iterator methodIter = methods.begin(); methodIter != methods.end(); ++methodIter)
	{
		tables.emplace(methodIter->first, std::move(Richardson()));
	}
}

void MethodWrapperBase::initalize(const OdeSolverParams& paramsIn)
{
	//Build the solvers
	buildSolvers(paramsIn);

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
		size_t currentVectorSize = findMethod(static_cast<SolverIF::SOLVER_TYPES>(richItr->first))->getCurrentState().size();

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

//Build up our solvers
void MethodWrapperBase::buildSolvers(const OdeSolverParams& paramsIn)
{

	//Check if our problem is stiff or large eigenvalues that might affect stability
	if (paramsIn.isStiff || paramsIn.isFast)
	{
		//Add Implict Methods only
		
		//Exit here
		return;
	}
	//Check problem is large computationally
	else if (paramsIn.isLarge)
	{
		//Add RK4 method
		getMethodMap().emplace(static_cast<unsigned int>(SolverIF::SOLVER_TYPES::RUNGE_KUTTA_FOUR),
			std::move(unique_ptr<SolverIF>(new RK4)));

		//Exit here
		return;
	}
	//What methods do we want to use
	else
	{
		//Add Euler Method
		if (paramsIn.useEuler)
		{
			//Add Euler to our allowed methods
			getMethodMap().emplace(static_cast<unsigned int>(SolverIF::SOLVER_TYPES::EULER),
				std::move(unique_ptr<SolverIF>(new Euler)));
		}

		//Add RK4
		if (paramsIn.useRK4)
		{
			//Add Runge Kutta 4 to our allowed methods
			getMethodMap().emplace(static_cast<unsigned int>(SolverIF::SOLVER_TYPES::RUNGE_KUTTA_FOUR),
				std::move(unique_ptr<SolverIF>(new RK4)));
		}
	}
}

void MethodWrapperBase::clearMethods()
{
	//Clear our methods
	methods.clear();

	//Clear our tables
	tables.clear();
}