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

void OdeSolver::runMethod(const OdeFunIF* problem, unique_ptr<SolverIF>& method, Richardson& tables, crvec initalCondition, rvec newState, const OdeSolverParams& currentParams, const double initalTime, const double newTime)
{
	//Loop over all the tables
	for (int i = 0; i < tables.getTableSize(); ++i)
	{
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

void OdeSolver::updateMethod(unique_ptr<SolverIF>& method, const OdeSolverParams& currentParams, Richardson& tables, crvec intialCondition, rvec newState, const double dt, const double initalTime, const double newTime, const OdeFunIF* problem, const int i)
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

void OdeSolver::buildSolution(unique_ptr<SolverIF>& currentMethod,const unsigned int currentMethodId, crvec initalCondition, const OdeFunIF* problem, const double beginTime, const double endTime)
{
	//Current method parameters
	OdeSolverParams& currentMethodParams = params.find(currentMethodId)->second;

	//Reset the satisfaction criteria
	currentMethodParams.satifiesError = false;

	//Get our richardson tables
	Richardson& currentTable = methods.get()->getTableMap().find(currentMethodId)->second;

	//Initalize our vector for the updated result
	vec newState;

	//Set our inital convergence criterial the the theoretical local truncation error 
	currentMethodParams.c = currentMethod->getErrorOrder() + static_cast<double>(currentMethodParams.minTableSize);

	//Set the currentError to a default value
	currentMethodParams.currentError = 99999.;

	//Update dt with our convergence criteria
	updateDt(currentMethodParams, true);

	//Get the current time
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

	//initalize our current starting table size to the min table size
	currentMethodParams.currentTableSize = generalParams.minTableSize;

	int counter = 0;

	//Run each result several times
	do
	{
		//Update our table
		currentTable.initalizeSteps(currentMethodParams.redutionFactor, currentMethodParams.dt);

		//Update the Richardson Table Size
		currentTable.BuildTables(currentMethodParams.currentTableSize, initalCondition.size());

		//Run our method
		runMethod(
			problem,
			currentMethod,
			currentTable,
			initalCondition,
			newState,
			currentMethodParams,
			beginTime,
			endTime);

		//Update the results with the new error
		currentMethodParams.currentError = currentTable.error(newState, currentMethodParams.c, currentMethod->getErrorOrder());

		//std::cout << currentMethodParams.c << "\n";

		//Build more tables if the error is greater then the greatest error
	} while (updateDt(currentMethodParams, false));// && currentMethodParams.currentTableSize++ < currentMethodParams.maxTableSize);

	//Get the second time point
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

	//Save the duriation of time
	currentMethodParams.currentRunTime = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count();
}

const bool OdeSolver::updateDt(OdeSolverParams& currentParams, const bool firstPassThrough)
{
	//Get our current error
	const double& currentError = currentParams.currentError;

	//Check to see if we are less than the greatest allowed error to show convergence
	if (currentError <= currentParams.upperError)
	{
		//We coverged to error we wanted
		currentParams.satifiesError = true;
		return false;
	}
	//Did not converge so update parameters
	else
	{
		//Get our dt
		double& dt = currentParams.dt;

		//Get our convergence estimate
		double& c = currentParams.c;

		//Estimate our current constant if not first pass through
		double funcDerv = 0.0;
		if (!firstPassThrough)
		{
			funcDerv = currentError / pow(dt, c);
		}
		else
		{
			funcDerv = 1.0;
		}

		//Set our clamped flag to false
		currentParams.isDtClamped = false;

		//Get our desired lower error
		const double& desiredErrorLeft = currentParams.lowerError;

		//Get our desired upper error
		const double& desiredErrorRight = currentParams.upperError;

		//Get our dt for lower desired error
		double leftDt = pow(desiredErrorLeft / funcDerv, 1. / c);

		//Get our dt for upper desired error
		double rightDt = pow(desiredErrorRight / funcDerv, 1. / c);

		//Get a dt for the middile
		dt = (rightDt + leftDt) / 2.;

		//Reset dt if it gets too small (the program would never finish)
		if (!isfinite(dt) && !firstPassThrough)
		{
			//Reset dt to the smallest parameters
			dt = currentParams.minDt;
		}
		//Dt is outside are min dt window we want to clamp it
		if (dt < currentParams.minDt && !firstPassThrough)
		{
			//Set the clamping flag
			currentParams.isDtClamped = true;

			//Update the table size to hope for better convergence
			currentParams.currentTableSize++;

			//Clamp dt
			dt = currentParams.minDt;
		}
		else if (dt > currentParams.maxDt && !firstPassThrough)
		{
			//Set the clamping flag
			currentParams.isDtClamped = true;

			//Clamp dt
			dt = currentParams.maxDt;

			//Update the table size to hope for better convergence
			currentParams.currentTableSize--;
		}
		else
		{
			//Nothing to do here
		}

		//We did not coverged to error we wanted
		currentParams.satifiesError = false;

		//If dt is too large we want to decrease the result
		if (currentParams.currentTableSize < currentParams.minTableSize)
		{
			//Clamp the table
			currentParams.currentTableSize = currentParams.minTableSize;
			return true;
		}
		//If dt is too small we dont want to keep growing
		else if (currentParams.currentTableSize > currentParams.maxTableSize)
		{
			//Reduce the table size if it is too large
			currentParams.currentTableSize--;
			return false;
		}
		else
		{
			return true;
		}
	}
}

void OdeSolver::run(OdeFunIF* problem, crvec initalConditions, const double beginTime, const double endTime)
{
	//Initalize all the methods
	methods.get()->updateAll(initalConditions, generalParams.minTableSize, generalParams.redutionFactor, generalParams.dt);

	//Get the method map
	methodMap& allowedMethods = methods->getMethodMap();

	//Clear out our threads from the previous run (will be added later)
	richardsonThreads.clear();

	//Iterate over all the methods
	for (methodMap::iterator methodItr = allowedMethods.begin(); methodItr != allowedMethods.end(); ++methodItr)
	{
		//Solver for the next time step for the current method
		//richardsonThreads.push_back(thread(&OdeSolver::buildSolution, methodItr->second, initalConditions, problem, beginTime, endTime));
		buildSolution(methodItr->second, methodItr->first, initalConditions, problem, beginTime, endTime);
	}

	//Join all the threads to get the results
	/*
	for (vector<thread>::iterator threadItr = richardsonThreads.begin(); threadItr != richardsonThreads.end(); ++threadItr)
	{
		threadItr->join();
	}
	*/

	std::cout << methods.get()->findMethod(SolverIF::SOLVER_TYPES::EULER)->getCurrentState()[0] << "\n";
	std::cout << methods.get()->findTable(SolverIF::SOLVER_TYPES::EULER).error() << "\n";
	std::cout << this->params.find(static_cast<unsigned int>(SolverIF::SOLVER_TYPES::EULER))->second.currentRunTime << "\n";
	std::cout << "\n\n" << this->params.find(static_cast<unsigned int>(SolverIF::SOLVER_TYPES::EULER))->second.dt << "\n";
	std::cout << "\n\n" << this->params.find(static_cast<unsigned int>(SolverIF::SOLVER_TYPES::EULER))->second.currentTableSize << "\n";
	std::cout << this->params.find(static_cast<unsigned int>(SolverIF::SOLVER_TYPES::EULER))->second.c << "\n";

	std::cout << methods.get()->findMethod(SolverIF::SOLVER_TYPES::RUNGE_KUTTA_FOUR)->getCurrentState()[0] << "\n";
	std::cout << methods.get()->findTable(SolverIF::SOLVER_TYPES::RUNGE_KUTTA_FOUR).error() << "\n";
	std::cout << this->params.find(static_cast<unsigned int>(SolverIF::SOLVER_TYPES::RUNGE_KUTTA_FOUR))->second.currentRunTime << "\n";
	std::cout << "\n\n" << this->params.find(static_cast<unsigned int>(SolverIF::SOLVER_TYPES::RUNGE_KUTTA_FOUR))->second.dt << "\n";
	std::cout << "\n\n" << this->params.find(static_cast<unsigned int>(SolverIF::SOLVER_TYPES::RUNGE_KUTTA_FOUR))->second.currentTableSize << "\n";
	std::cout << this->params.find(static_cast<unsigned int>(SolverIF::SOLVER_TYPES::RUNGE_KUTTA_FOUR))->second.c << "\n";
}