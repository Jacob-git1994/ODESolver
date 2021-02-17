#include "MethodWrapperBase.h"

/// <summary>
/// Finds the method by the methods enumeration.
/// </summary>
/// <param name="solver"></param>
/// <returns></returns>
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

/// <summary>
/// Finds the table based on the method enumeration.
/// </summary>
/// <param name="solver"></param>
/// <returns></returns>
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

/// <summary>
/// Initalizes all the methods with the size of our vector.
/// </summary>
/// <param name="state"></param>
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

/// <summary>
/// Get the method map.
/// </summary>
/// <returns></returns>
methodMap& MethodWrapperBase::getMethodMap()
{
	return methods;
}

/// <summary>
/// Get the table map.
/// </summary>
/// <returns></returns>
tableMap& MethodWrapperBase::getTableMap()
{
	return tables;
}

/// <summary>
/// Get a const method map.
/// </summary>
/// <returns></returns>
const methodMap& MethodWrapperBase::getMethodMap() const
{
	return methods;
}

/// <summary>
/// Get a const table map.
/// </summary>
/// <returns></returns>
const tableMap& MethodWrapperBase::getTableMap() const
{
	return tables;
}

/// <summary>
/// Build up our richardson tables for each allowed method
/// </summary>
void MethodWrapperBase::buildTables()
{
	//Go through what methods we are using and add a table class to each one
	for (methodMap::const_iterator methodIter = methods.cbegin(); methodIter != methods.cend(); ++methodIter)
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

	//Check if our problem is stiff
	if (paramsIn.isStiff)
	{
		//Add Implict Methods only
		
		//Exit here
		return;
	}
	//Check problem might grow really fast use rk4 & implict methods
	else if (paramsIn.isFast)
	{
		//Add RK4 method
		getMethodMap().emplace(static_cast<unsigned int>(SolverIF::SOLVER_TYPES::RUNGE_KUTTA_FOUR),
			std::move(unique_ptr<SolverIF>(new RK4)));
	}
	//Check if the system is large so the run time could be large with implict but we need high accuracy
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

		//Add RK2
		if (paramsIn.useRK2)
		{
			getMethodMap().emplace(static_cast<unsigned int>(SolverIF::SOLVER_TYPES::RUNGE_KUTTA_TWO),
				std::move(unique_ptr<SolverIF>(new RK2)));
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