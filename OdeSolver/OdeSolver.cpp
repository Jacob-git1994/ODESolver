#include "OdeSolver.h"

OdeSolver::OdeSolver(const OdeSolverParams& paramsIn) :
	generalParams(paramsIn)
{
	setup();
}

void OdeSolver::runMethod(const OdeFunIF* problem, unique_ptr<SolverIF>& method, Richardson& tables, crvec initalCondition, rvec newState, const OdeSolverParams& currentParams, const double initalTime, const double newTime)
{
	//Loop over all the tables
	for (unsigned int i = 0; i < tables.getTableSize(); ++i)
	{
		//Update the tables
		updateMethod(method, currentParams, tables, initalCondition, newState, currentParams.dt, initalTime, newTime, problem, i);
	}
}

void OdeSolver::updateMethod(unique_ptr<SolverIF>& method, const OdeSolverParams& currentParams, Richardson& tables, crvec intialCondition, rvec newState, const double dt, const double initalTime, const double newTime, const OdeFunIF* problem, const int i)
{
	//Solve for the next time step
	method->update(intialCondition, newState, currentParams.dt / pow(tables.getReductionFactor(), static_cast<double>(i)), initalTime, static_cast<int>(pow(tables.getReductionFactor(), static_cast<double>(i))), problem);

	//Add the result into the tables
	tables(i, 0, newState);
}

/// <summary>
/// Runs one step for our time stepping algorithem. We iterativly redo each step measuring the convergence and determine a new eastimate for dt
/// in order to keep the error bounded
/// </summary>
/// <param name="currentMethod"></param>
/// <param name="currentMethodId"></param>
/// <param name="initalCondition"></param>
/// <param name="problem"></param>
/// <param name="beginTime"></param>
/// <param name="endTime"></param>
/// <returns></returns>
vec OdeSolver::buildSolution(unique_ptr<SolverIF>& currentMethod,const unsigned int currentMethodId, Richardson& currentTable, OdeSolverParams& currentMethodParams, crvec initalCondition, const OdeFunIF* problem, const double beginTime, const double endTime)
{
	//Reset the satisfaction criteria
	currentMethodParams.satifiesError = false;

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
	} while (updateDt(currentMethodParams, false, beginTime, endTime));

	//Get the second time point
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

	//Save the duriation of time
	currentMethodParams.currentRunTime += std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count();

	//Return the new state found
	return newState;
}

/// <summary>
/// We check if dt is valid and supports the error goal.
/// If dt fails the checks we find a new dt based on the convergence estimate. If we are on the last time step, we clamp dt so the
/// last step run will be at the final end time
/// </summary>
/// <param name="currentParams"></param>
/// <param name="firstPassThrough"></param>
/// <param name="beginTime"></param>
/// <param name="endTime"></param>
/// <returns></returns>
const bool OdeSolver::updateDt(OdeSolverParams& currentParams, const bool firstPassThrough, const double beginTime, const double endTime)
{
	//Get our variables to use to upgrade dt
	double& dt = currentParams.dt;
	const double& covergenceEstimate = currentParams.c;
	const double& currentError = currentParams.currentError;
	const double& desiredError = currentParams.upperError;
	const double& lowestErrorAllowed = currentParams.lowerError;
	const double& minDtUpgrade = currentParams.minDt;
	const double& maxDtUpgrade = currentParams.maxDt;
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
	}

	//If this is first pass through we dont want to break
	if (firstPassThrough && !currentParams.lastRun)
	{
		//Set the current table size
		currentParams.currentTableSize = currentParams.minTableSize;
	}
	//We have converged within desired error range
	else if ((!firstPassThrough && currentError <= desiredError && currentError >= lowestErrorAllowed && !currentParams.lastRun) || (!isfinite(covergenceEstimate) || covergenceEstimate < 0.0))
	{
		//Break our loop
		return false;
	}
	//We did not converge to desired solution
	else if (!currentParams.lastRun)
	{
		//Set a guess for the derivative
		funcDerv = currentError / pow(dt, covergenceEstimate);

		//Get an estimate on how much we want to increase/decrease dt
		desiredUpgrade = pow(desiredError / (funcDerv * currentError), 1. / covergenceEstimate);
		desiredUpgrade = std::max(desiredUpgrade, minDtUpgrade);
		desiredUpgrade = std::min(desiredUpgrade, maxDtUpgrade);

		//Update our dt
		dt *= .9 * desiredUpgrade;

		//Check how we want to modify our table size
		if (currentError < lowestErrorAllowed)
		{
			currentParams.currentTableSize--;
		}
		else if (currentError > desiredError)
		{
			currentParams.currentTableSize++;
		}
		else
		{
			//Do nothing here
		}

		//Check to ensure we clamp table size
		if (currentParams.currentTableSize > currentParams.maxTableSize)
		{
			currentParams.currentTableSize = currentParams.maxTableSize;
		}
		else if (currentParams.currentTableSize < currentParams.minTableSize)
		{
			currentParams.currentTableSize = currentParams.minTableSize;
		}
		else
		{
			//Do nothing here
		}
	}
	//We ran the last time
	else
	{
		return false;
	}

	return true;
}

/// <summary>
/// This runs the core algorithem for each allowed method
/// </summary>
/// <param name="problem"></param>
/// <param name="initalConditions"></param>
/// <param name="beginTime"></param>
/// <param name="endTime"></param>
/// <param name="nodes"></param>
void OdeSolver::run(OdeFunIF* problem, crvec initalConditions, const double beginTime, const double endTime, const unsigned int nodes)
{
	//Initalize all the methods
	methods.updateAll(initalConditions, generalParams.minTableSize, generalParams.redutionFactor, generalParams.dt);

	//Get the method map
	methodMap& allowedMethods = methods.getMethodMap();

	//Check if we have avaiable methods
	if (allowedMethods.empty())
	{
		throw runtime_error("No Allowed Methods Available");
	}
	else
	{
		//Clear out our threads from the previous run (will be added later)
		methodThreads.clear();

		//Iterate over all the methods
		for (methodMap::iterator methodItr = allowedMethods.begin(); methodItr != allowedMethods.end(); ++methodItr)
		{
			//Get the current method
			unique_ptr<SolverIF>& currentMethod = methodItr->second;

			//Get the current parameters associated with the method
			OdeSolverParams& currentParams = params.find(methodItr->first)->second;

			//Get the current richardson parameters
			Richardson& currentTables = methods.getTableMap().find(methodItr->first)->second;

			//Put back the time iterations
			methodThreads.push_back(std::move(thread(
				&OdeSolver::updateNextTimeStep,
				this, 
				methodItr->first, 
				std::ref(currentMethod),
				std::ref(currentParams),
				std::ref(currentTables),
				beginTime, 
				endTime, 
				initalConditions, 
				std::cref(problem))));
		}

		//Join all the threads to get the results
		for (vector<thread>::iterator threadItr = methodThreads.begin(); threadItr != methodThreads.end(); ++threadItr)
		{
			threadItr->join();
		}
	}

	std::cout << setprecision(15) << methods.findMethod(SolverIF::SOLVER_TYPES::EULER)->getCurrentState()[0] << "\n";
	std::cout << methods.findTable(SolverIF::SOLVER_TYPES::EULER).error() << "\n";
	std::cout << this->params.find(static_cast<unsigned int>(SolverIF::SOLVER_TYPES::EULER))->second.currentRunTime << "\n";
	std::cout << "\n\n" << this->params.find(static_cast<unsigned int>(SolverIF::SOLVER_TYPES::EULER))->second.dt << "\n";
	std::cout << "\n\n" << this->params.find(static_cast<unsigned int>(SolverIF::SOLVER_TYPES::EULER))->second.currentTableSize << "\n";
	std::cout << this->params.find(static_cast<unsigned int>(SolverIF::SOLVER_TYPES::EULER))->second.c << "\n";

	
	std::cout << setprecision(15) << methods.findMethod(SolverIF::SOLVER_TYPES::RUNGE_KUTTA_FOUR)->getCurrentState()[0] << "\n";
	std::cout << methods.findTable(SolverIF::SOLVER_TYPES::RUNGE_KUTTA_FOUR).error() << "\n";
	std::cout << this->params.find(static_cast<unsigned int>(SolverIF::SOLVER_TYPES::RUNGE_KUTTA_FOUR))->second.currentRunTime << "\n";
	std::cout << "\n\n" << this->params.find(static_cast<unsigned int>(SolverIF::SOLVER_TYPES::RUNGE_KUTTA_FOUR))->second.dt << "\n";
	std::cout << "\n\n" << this->params.find(static_cast<unsigned int>(SolverIF::SOLVER_TYPES::RUNGE_KUTTA_FOUR))->second.currentTableSize << "\n";
	std::cout << this->params.find(static_cast<unsigned int>(SolverIF::SOLVER_TYPES::RUNGE_KUTTA_FOUR))->second.c << "\n";
	
}

/// <summary>
/// When the object is already initalized and we want to update the parameters for another run we clear out our maps
/// and rebuild the methods with the new parameters
/// </summary>
/// <param name="paramsIn"></param>
void OdeSolver::refreshParams(const OdeSolverParams& paramsIn)
{
	//Cleaar our methods
	methods.clearMethods();

	//Clear out our method parameters
	params.clear();

	//Clear out our results
	resultMap.clear();

	//Clear out our vector of threads
	methodThreads.clear();

	//ReInitaize our general parameters
	generalParams = paramsIn;

	//Set up everything again
	setup();
}

/// <summary>
/// Build up all our maps for our methods
/// </summary>
void OdeSolver::setup()
{
	//Check the user inputs
	if (!generalParams.checkUserInputs())
	{
		throw invalid_argument("Invalid Ode Parameters");
		exit(-1);
	}

	//Build our allowed methods based on the parameter input
	try
	{
		methods.initalize(generalParams);
	}
	catch (exception& e)
	{
		cerr << e.what();
		exit(-1);
	}

	//Get the methods allowed
	const methodMap& allowedMethods = methods.getMethodMap();

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
		params.emplace(methodId, generalParams);
	}
}

/// <summary>
/// Build the current methods solution until the final time specified.
/// </summary>
/// <param name="methodId"></param>
/// <param name="currentMethod"></param>
/// <param name="currentParameters"></param>
/// <param name="currentTables"></param>
/// <param name="beginTime"></param>
/// <param name="endTime"></param>
/// <param name="initalConditions"></param>
/// <param name="problem"></param>
void OdeSolver::updateNextTimeStep(
	const unsigned int methodId, 
	unique_ptr<SolverIF>& currentMethod, 
	OdeSolverParams& currentParameters, 
	Richardson& currentTables, 
	const double beginTime, 
	const double endTime, 
	const valarray<double>& initalConditions, 
	const OdeFunIF* problem)
{
	//Get this current methods dt
	double& dt = currentParameters.dt;

	//Save our current time
	double currentTime = beginTime;

	//Save our currentState
	valarray<double> currentState = initalConditions;

	while (currentTime < endTime)
	{
		//Solver for the next time step for the current method
		currentState = buildSolution(currentMethod, methodId, currentTables, currentParameters, currentState, problem, currentTime, endTime);

		//Update the time
		currentTime += dt;

		//In future I plan on adding each result to our solution map
	}
}