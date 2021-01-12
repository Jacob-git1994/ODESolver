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

void OdeSolver::runMethod(OdeFunIF* problem, unique_ptr<SolverIF>& method, Richardson& tables, crvec initalCondition, rvec newState, const OdeSolverParams& currentParams, const double initalTime, const double newTime)
{
	//Loop over all the tables
	for (int i = 0; i < tables.getTableSize(); ++i)
	{
		//Solve for the next time step in parallel
		/*Does not support this yet..
		richardsonThreads.push_back(thread(&OdeSolver::updateMethod, 
			std::ref(method), 
			currentParams, 
			std::ref(tables), 
			initalCondition, 
			std::ref(newState), 
			currentParams.dt, 
			initalTime, newTime, 
			std::ref(problem), 
			i));
		*/

		//Update the tables
		updateMethod(method, currentParams, tables, initalCondition, newState, currentParams.dt, initalTime, newTime, problem, i);
	}

	//Join everything back together (Not implimented yet)
	/*
	for (threadVector::iterator vectorItr = richardsonThreads.begin(); vectorItr != richardsonThreads.end(); ++vectorItr)
	{
		vectorItr->join();
	}
	*/
}

void OdeSolver::updateMethod(unique_ptr<SolverIF>& method, const OdeSolverParams& currentParams, Richardson& tables, crvec intialCondition, rvec newState, const double dt, const double initalTime, const double newTime, OdeFunIF* problem, const int i)
{
	//Solve for the next time step
	method->update(intialCondition, newState, currentParams.dt / pow(tables.getReductionFactor(), static_cast<double>(i)), initalTime, newTime, problem);

	//Add the result into the tables
	tables(i, 0, newState);
}

void OdeSolver::gatherParameters(OdeSolverParams& currentParams, Richardson& currentTable, const unsigned int& currentMethod)
{
	//Get our current parameters
	currentParams = params.find(currentMethod)->second;

	//Get our current table
	currentTable = methods.get()->getTableMap().find(currentMethod)->second;
}

void OdeSolver::buildSolution(crvec initalCondition, OdeFunIF* problem, const double beginTime, const double endTime)
{
	//Initalize all the methods
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

		//Reset the satisfaction criteria
		currentMethodParams.satifiesError = false;

		//Get our richardson tables
		Richardson& currentTable = methods.get()->getTableMap().find(currentMethodId)->second;

		//Initalize our vector for the updated result
		vec newState;

		//Clear out our threads from the previous run (will be added later)
		richardsonThreads.clear();

		//Set dt to the max dt allowed
		currentMethodParams.dt = currentMethodParams.maxDt;

		//Get the current time
		std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

		//initalize our current starting table size to the min table size
		currentMethodParams.currentTableSize = generalParams.minTableSize;

		//Update our table
		currentTable.initalizeSteps(currentMethodParams.redutionFactor, currentMethodParams.dt);

		//Run each result several times
		do
		{
			//Update the Richardson Table Size
			currentTable.BuildTables(currentMethodParams.currentTableSize, initalCondition.size());

			//Run our method
			runMethod(
				problem,
				methodItr->second,
				currentTable,
				initalCondition,
				newState,
				currentMethodParams,
				beginTime,
				endTime);

			//Update the results with the new error
			currentMethodParams.currentError = currentTable.error(newState, currentMethodParams.c, methodItr->second->getErrorOrder());

			//std::cout << currentMethodParams.c << "\n";

			//Build more tables if the error is greater then the greatest error
		} while (updateDt(currentMethodParams));

		//Correct the current table size of its last increment
		currentMethodParams.currentTableSize--;

		//Get the second time point
		std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

		//Save the duriation of time
		currentMethodParams.currentRunTime = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count();
	}
}

const bool OdeSolver::updateDt(OdeSolverParams& currentParams)
{
	//Get our current error
	const double& currentError = currentParams.currentError;

	//Check to see if our error has been satisfied.
	if (currentError >= currentParams.lowerError && currentError <= currentParams.upperError)
	{
		return false;
	}
	//Did not converge so update parameters
	else
	{
		//Get our dt
		double& dt = currentParams.dt;

		//Get our convergence estimate
		double& c = currentParams.c;

		//Get our desired lower error
		const double& desiredErrorLeft = currentParams.lowerError;

		//Get our desired upper error
		const double& desiredErrorRight = currentParams.upperError;

		//Get our dt for lower desired error
		double leftDt = pow(desiredErrorLeft, 1. / c);

		//Get our dt for upper desired error
		double rightDt = pow(desiredErrorRight, 1. / c);

		//Get a dt for the middile
		dt = (rightDt - leftDt) / 2.;

		return true;
	}
}

void OdeSolver::run(OdeFunIF* problem, crvec initalConditions, const double beginTime, const double endTime)
{
	//Solver for the next time step
	buildSolution(initalConditions, problem, beginTime, endTime);

	std::cout << methods.get()->findMethod(SolverIF::SOLVER_TYPES::EULER)->getCurrentState()[0] << "\n";
	std::cout << methods.get()->findTable(SolverIF::SOLVER_TYPES::EULER).error() << "\n";
	std::cout << this->params.find(static_cast<unsigned int>(SolverIF::SOLVER_TYPES::EULER))->second.currentRunTime << "\n";
	std::cout << "\n\n" << this->params.find(static_cast<unsigned int>(SolverIF::SOLVER_TYPES::EULER))->second.dt;
	std::cout << "\n\n" << this->params.find(static_cast<unsigned int>(SolverIF::SOLVER_TYPES::EULER))->second.currentTableSize;
	std::cout << this->params.find(static_cast<unsigned int>(SolverIF::SOLVER_TYPES::EULER))->second.c;
}