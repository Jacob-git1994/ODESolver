#include "OdeSolver.h"

bool OdeSolver::isAllExplict() const
{
	return generalParams.useEuler && generalParams.useRK2 && generalParams.useRK4;
}

bool OdeSolver::isAllImplict() const
{
	return generalParams.useCrank && generalParams.useImplictEuler;
}

OdeSolver::OdeSolver(const OdeSolverParams& paramsIn) :
	generalParams(paramsIn),
	methods(nullptr)
{
	//Check the user inputs
	if (!generalParams.checkUserInputs())
	{
		throw invalid_argument("Invalid Ode Parameters");
		exit(-1);
	}

	//This section will build up what methods we want to attempt to use
	//Case if all methods are explictly defined
	if (isAllExplict() && !isAllImplict()) 
	{
		//We shouldn't get here...yet
	}
	//Case if all methiods are explictly defined
	else if (isAllImplict() && !isAllExplict()) 
	{
		//we shouldn't get here...yet
	}
	else
	{
		//Use only explict methods
		if (generalParams.isLarge || generalParams.isFast)
		{
			//we shouldnt get here yet
		}
		//Use only implict methods
		else if (generalParams.isStiff)
		{

		}
		//Case where we just want to build up all the fun methods
		else
		{
			methods = unique_ptr<MethodWrapperIF>(new MethodWrapper);
		}
	}

	//Check to see if methods was even allocated correctly
	if (!methods)
	{
		throw bad_alloc();
	}

	//Initalize our methods
	try
	{
		methods.get()->initalize();
	}
	catch (exception& e)
	{
		cerr << e.what();
		exit(-1);
	}

	//Get the methods allowed
	const methodMap& allowedMethods = methods.get()->getMethodMap();

	//Check if no methods were created
	if (allowedMethods.empty())
	{
		throw runtime_error("Method map is zero");
		exit(-1);
	}

	//Iterate through the allowed methods to build the parameter tables
	for (methodMap::const_iterator methodItr = allowedMethods.cbegin(); methodItr != allowedMethods.cend(); ++methodItr)
	{
		//This is the ID/enum for the method
		const unsigned int methodId = methodItr->first;

		//Add the method id and parameters to our param map
		params.emplace(methodId, paramsIn);
	}
}

void OdeSolver::findOptimalParams(crvec initalCondition,OdeFunIF* problem, const double currentTime, const double newTime)
{
	//Initalize our method wrapper
	methods.get()->updateAll(initalCondition, generalParams.minTableSize, generalParams.redutionFactor, generalParams.dt);

	//Get the method map
	methodMap& allowedMethods = methods.get()->getMethodMap();

	//Iterate over all the methods
	for (methodMap::iterator methodItr = allowedMethods.begin(); methodItr != allowedMethods.end(); ++methodItr)
	{

		//Current method id
		unsigned int currentMethodId = methodItr->first;

		//Current method parameters
		OdeSolverParams& currentMethodParams = params.find(currentMethodId)->second;

		//Get our richardson tables
		Richardson& currentTable = methods.get()->getTableMap().find(currentMethodId)->second;

		//Set our table size
		currentMethodParams.currentTableSize = generalParams.minTableSize;

		//Initalize our vector for the updated result
		vec newState;

		//Set our current run time to the min run time
		currentMethodParams.currentRunTime = currentMethodParams.minRunTime;

		//Run until we break the max run time
		do
		{
			//Get the current time
			std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

			//Run through building all tables
			do
			{
				//Update the Richardson Table Size
				currentTable.BuildTables(currentMethodParams.currentTableSize, initalCondition.size());

				//Run our method
				runMethod(
					problem,
					methodItr->second.get(),
					currentTable,
					initalCondition,
					newState,
					currentMethodParams,
					currentTime,
					newTime);

				//Update the results with the new error
				currentMethodParams.currentError = currentTable.error();

			} while (currentMethodParams.currentTableSize++ <= generalParams.maxTableSize &&
				(currentMethodParams.currentError <= currentMethodParams.lowerError ||
					currentMethodParams.currentError >= currentMethodParams.upperError));

			//Get the second time point
			std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

			//Save the duriation of time
			currentMethodParams.currentRunTime = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count();

			//Update dt
			currentMethodParams.dt /= currentMethodParams.redutionFactor;

		} while (currentMethodParams.currentRunTime < currentMethodParams.maxRunTime && 
			(currentMethodParams.currentError <= currentMethodParams.lowerError ||
			currentMethodParams.currentError >= currentMethodParams.upperError));
	}
}

void OdeSolver::runMethod(OdeFunIF* problem, SolverIF* method, Richardson& tables, crvec initalCondition, rvec newState, const OdeSolverParams& currentParams, const double initalTime, const double newTime)
{
	//Loop over all the tables
	for (int i = 0; i < tables.getTableSize(); ++i)
	{
		//Solve for the next time step
		method->update(initalCondition,newState, currentParams.dt/pow(tables.getReductionFactor(),static_cast<double>(i)), initalTime, newTime, problem);

		//Add the result into the tables
		tables(i, 0, newState);
	}
}

void OdeSolver::run(OdeFunIF* problem,crvec initalConditions,const double beginTime, const double endTime)
{
	//Find the best set of parameters (This will be optimized later
	findOptimalParams(initalConditions, problem, beginTime, endTime);

	std::cout << methods.get()->findMethod(SolverIF::SOLVER_TYPES::EULER)->getCurrentState()[0] << "\n";
	std::cout << methods.get()->findTable(SolverIF::SOLVER_TYPES::EULER).error() << "\n";
	std::cout << this->params.find(static_cast<unsigned int>(SolverIF::SOLVER_TYPES::EULER))->second.currentRunTime;
}