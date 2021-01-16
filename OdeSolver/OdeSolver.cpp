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
	method->update(intialCondition, newState, currentParams.dt / pow(tables.getReductionFactor(), static_cast<double>(i)), initalTime, static_cast<int>(pow(tables.getReductionFactor(), static_cast<double>(i))), problem);

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

vec OdeSolver::buildSolution(unique_ptr<SolverIF>& currentMethod,const unsigned int currentMethodId, crvec initalCondition, const OdeFunIF* problem, const double beginTime, const double endTime)
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

	//Update dt with our convergence criteria
	updateDt(currentMethodParams, true, beginTime, endTime);

	//Get the current time
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

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

		//Build more tables if the error is greater then the greatest error
	} while (updateDt(currentMethodParams, false, beginTime, endTime));// && currentMethodParams.currentTableSize++ < currentMethodParams.maxTableSize);

	//Get the second time point
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

	//Save the duriation of time
	currentMethodParams.currentRunTime = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count();

	//Return the new state found
	return newState;
}

const bool OdeSolver::updateDt(OdeSolverParams& currentParams, const bool firstPassThrough, const double beginTime, const double endTime)
{
	//Get our variables to use to upgrade dt
	double& dt = currentParams.dt;
	const double& covergenceEstimate = currentParams.c;
	const double& currentError = currentParams.currentError;
	const double& desiredError = currentParams.upperError;
	const double& minDtUpgrade = currentParams.minDt;
	const double& maxDtUpgrade = currentParams.minDt;
	double funcDerv = 0.0;
	double desiredUpgrade = 0.0;

	//Reset dt if it would fall outside our bounds
	if (dt + beginTime >= endTime && !currentParams.lastRun)
	{
		//Set dt to the remaining difference
		dt = endTime - beginTime;

		//Set the last run flag
		currentParams.lastRun = true;

		//Set the table size to max to try and get a good convergence
		currentParams.currentTableSize = currentParams.maxTableSize;

		//Run one more iteration
		return true;
	}

	//If this is first pass through we dont want to break
	if (firstPassThrough && !currentParams.lastRun)
	{
		//Set the current table size
		currentParams.currentTableSize = currentParams.minTableSize;
	}
	//We have converged
	else if (!firstPassThrough && currentError <= desiredError && !currentParams.lastRun)
	{
		//Break our loop
		return false;
	}
	//We did not converge to desired solution
	else if (!currentParams.lastRun)
	{
		//Set a guess for the derivative
		funcDerv = currentError / pow(dt, covergenceEstimate);

		//Get an estimate on how much we want to increase dt
		desiredUpgrade = pow(desiredError / (funcDerv * currentError), 1. / covergenceEstimate);
		desiredUpgrade = std::max(desiredUpgrade, 1.0 / currentParams.redutionFactor);
		desiredUpgrade = std::min(desiredUpgrade, maxDtUpgrade);

		//Update our dt
		dt *= .9 * desiredUpgrade;

		//Update the table size
		currentParams.currentTableSize++;

		//Check to ensure we clamp table size
		if (currentParams.currentTableSize > currentParams.maxTableSize)
		{
			currentParams.currentTableSize = currentParams.maxTableSize;
		}
		else if (currentParams.currentTableSize < currentParams.minTableSize)
		{
			currentParams.currentTableSize = currentParams.minTableSize;
		}

		return true;
	}
	//We ran the last time
	else
	{
		return false;
	}

	return true;
}

void OdeSolver::run(OdeFunIF* problem, crvec initalConditions, const double beginTime, const double endTime, const unsigned int nodes)
{
	//Initalize all the methods
	methods.get()->updateAll(initalConditions, generalParams.minTableSize, generalParams.redutionFactor, generalParams.dt);

	//Get the method map
	methodMap& allowedMethods = methods->getMethodMap();

	//Clear out our threads from the previous run (will be added later)
	methodThreads.clear();

	//Iterate over all the methods
	for (methodMap::iterator methodItr = allowedMethods.begin(); methodItr != allowedMethods.end(); ++methodItr)
	{
		//Save our current time
		double currentTime = beginTime;

		//Save our currentState
		valarray<double> currentState = initalConditions;

		//Get this methods dt
		double& dt = params.find(methodItr->first)->second.dt;

		while (currentTime < endTime)
		{
			//Solver for the next time step for the current method
			currentState = buildSolution(methodItr->second, methodItr->first, currentState, problem, currentTime, endTime);

			//Update the time
			currentTime += dt;
		}
	}

	std::cout << setprecision(15) << methods.get()->findMethod(SolverIF::SOLVER_TYPES::EULER)->getCurrentState()[0] << "\n";
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