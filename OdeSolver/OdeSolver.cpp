#include "OdeSolver.h"

/// <summary>
/// This constructor calls the set up method which builds up all our tables and maps to later be used when we decide to run.
/// </summary>
/// <param name="paramsIn"></param>
OdeSolver::OdeSolver(const OdeSolverParams& paramsIn) :
	generalParams(paramsIn)
{
	setup();
}

/// <summary>
/// Here we build up the richardson tables by calling the method solver i number of times.
/// Each ith run we divide up the stepsize by that much. 
/// Update method will then append the result into the richardson tables after a result is found. 
/// </summary>
/// <param name="problem"></param>
/// <param name="method"></param>
/// <param name="tables"></param>
/// <param name="initalCondition"></param>
/// <param name="newState"></param>
/// <param name="currentParams"></param>
/// <param name="initalTime"></param>
/// <param name="newTime"></param>
void OdeSolver::runMethod(const OdeFunIF* problem, unique_ptr<SolverIF>& method, const unsigned int currentMethodId, Richardson& tables, crvec initalCondition, rvec newState, const OdeSolverParams& currentParams, const double initalTime, const double newTime)
{
	//Loop over all the tables
	for (unsigned int i = 0; i < tables.getTableSize(); ++i)
	{
		//Check if we are using an implict method
		if (isExplict(currentMethodId))
		{
			//Solve for the next time step
			method->update(initalCondition, newState, currentParams.dt / pow(tables.getReductionFactor(), static_cast<double>(i)), initalTime, static_cast<int>(pow(tables.getReductionFactor(), static_cast<double>(i))), problem);
		}
		//If we are implict then run the implict updating method
		else if (!isExplict(currentMethodId))
		{
			//Solve for the next time step
			method->update(initalCondition, newState, currentParams.dt / pow(tables.getReductionFactor(), static_cast<double>(i)), initalTime, static_cast<int>(pow(tables.getReductionFactor(), static_cast<double>(i))), problem, currentParams.implictDt, currentParams.implictError);
		}
		else
		{
			throw runtime_error("Could not resolve method");
		}

		//Add the result into the tables
		tables.append(i, 0, newState);
	}
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
vec OdeSolver::buildSolution(unique_ptr<SolverIF>& currentMethod, const unsigned int currentMethodId, Richardson& currentTable, OdeSolverParams& currentMethodParams, crvec initalCondition, const OdeFunIF* problem,const double beginTime, const double endTime)
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
		runMethod(problem, currentMethod, currentMethodId, currentTable, initalCondition, newState, currentMethodParams, beginTime, endTime);

		//Update the results with the new error
		currentMethodParams.currentError = currentTable.error(newState, currentMethodParams.c);

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

	//Get the parameters we will modify
	double& dt = currentParams.dt;
	double& totalError = currentParams.totalError;
	double& upgradeFactor = currentParams.upgradeFactor;
	bool& clamp = currentParams.isDtClamped;
	bool& lastRun = currentParams.lastRun;
	bool& conditionsSatisfied = currentParams.satifiesError;
	size_t& currentTableSize = currentParams.currentTableSize;

	//Get parameters we will use
	const double& c = currentParams.c;
	const double& currentError = currentParams.currentError;
	const double& desiredError = currentParams.upperError;
	const double& lowestAllowableError = currentParams.lowerError;
	const double& minDtUpgrade = currentParams.minDt;
	const double& maxDtUpgrade = currentParams.maxDt;
	const double& smallestDtAllowed = currentParams.smallestAllowableDt;
	const size_t& maxTableSize = currentParams.maxTableSize;
	const size_t& minTableSize = currentParams.minTableSize;
	const bool& isStiff = currentParams.isStiff;
	const bool& isFast = currentParams.isFast;

	//Our desired upgrade ammount
	double desiredUpdate = 0.0;

	//Estimate the global errors (known error on grid + potential error left)
	const double globalError = totalError + std::floor((endTime - beginTime) / dt) * currentError;

	//If this is the first pass through want to reset our working parameters and set up our new step sizes
	if (firstPassThrough && !lastRun)
	{
		//Check if we want to start our search from the begining table size or start where we left off to increase performance
		if (!(isStiff || isFast || clamp))
		{
			currentTableSize = minTableSize;
		}

		//Update our dt by the upgrade factor found in previous iteration. If a previous iteration does not exsit then we skip this processing.
		if (upgradeFactor > 1.0)
		{
			dt *= upgradeFactor;
		}

		//Reset our parameters
		conditionsSatisfied = false;
		lastRun = false;
		clamp = false;

		//Check if we need to clamp dt if we are at the end point of the interval
		if (dt + beginTime > endTime)
		{
			//Reset dt to end where we plan on it ending
			dt = endTime - beginTime;

			//Update table size to max to hope for better convergence since we don't control dt anymore
			currentTableSize = maxTableSize;

			//Update the last run flag
			lastRun = true;
		}

		//Exit further processing
		return true;
	}
	else if (!lastRun)
	{
		//Check if we did not satisify the error and dt is not clampped
		if (globalError > desiredError && !lastRun && isfinite(c) && c > 0.0 && !clamp)
		{
			//Get an estimate on how much we want to increase/decrease dt
			desiredUpdate = pow(desiredError / globalError, 1. / c);
			desiredUpdate = std::max(desiredUpdate, minDtUpgrade);
			desiredUpdate = std::min(desiredUpdate, maxDtUpgrade);

			//Update the dt by the ratio we want
			dt *= .9 * desiredUpdate;

			//Increase our table size to increase accuracy
			currentTableSize++;
			
			//Clamp our table size if we exceed the bounds
			if (currentTableSize > maxTableSize)
			{
				//Reset the current table size
				currentTableSize = maxTableSize;
			}

			//Check to make sure dt satisfies what we will allow
			if (dt < smallestDtAllowed)
			{
				//Reset dt
				dt = smallestDtAllowed;

				//Update our flag
				clamp = true;
			}

			//Set our conditions satisifed to false as we did not converge to the correct solution
			conditionsSatisfied = false;

			return true;
		}
		//Check if we satisifed the error or dt was forced to be clampped we want to exit the iteration
		else if ((globalError <= desiredError || !isfinite(c) || c < 0.0 || clamp) && !lastRun)
		{
			//Get an estimate on how much we want to increase/decrease dt
			desiredUpdate = pow(desiredError / (currentError), 1. / c);
			desiredUpdate = std::max(desiredUpdate, minDtUpgrade);
			desiredUpdate = std::min(desiredUpdate, maxDtUpgrade);

			//Set our upgrade factor for the next run
			upgradeFactor = desiredUpdate;

			//Check if we can lower our table size
			if (globalError <= lowestAllowableError)
			{
				//Get an estimate on how much we want to increase/decrease dt base on if the error is too small
				desiredUpdate = pow(lowestAllowableError / (currentError), 1. / c);
				desiredUpdate = std::max(desiredUpdate, minDtUpgrade);
				desiredUpdate = std::min(desiredUpdate, maxDtUpgrade);

				//Change our upgrade factor
				upgradeFactor = desiredUpdate;

				//Decrease our table size
				currentTableSize--;
			}

			//Clamp our table size
			if (currentTableSize < minTableSize)
			{
				currentTableSize = minTableSize;
			}

			//Set our conditions satisfied
			conditionsSatisfied = true;

			//Accumulate the error
			currentParams.totalError += currentError;
			
			//Exit Processing and move onto the next iteration in time
			return false;
		}
	}
	//We ran with the last runs dt and now we can stop processing this method
	else
	{
		//check if we satisfied some error
		conditionsSatisfied = currentError <= desiredError;

		//Accumulate the error
		currentParams.totalError += currentError;

		//Exit processing and evaluate final result
		return false;
	}
	//If something goes wrong and we get to this step we want to just do it again :)
	return true;
}

/// <summary>
/// This runs the core alogirthm for all methods.
/// We push back each method to run in parallel on each thread. 
/// </summary>
/// <param name="problem"></param>
/// <param name="initalConditions"></param>
/// <param name="beginTime"></param>
/// <param name="endTime"></param>
/// <param name="nodes"></param>
void OdeSolver::run(const OdeFunIF* problem, crvec initalConditions, const double beginTime, const double endTime)
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
		//Clear out our threads from the previous run
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

			//Get the current problems result map
			vector<StateVector>& currentStateVector = resultMap.find(methodItr->first)->second;

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
				std::cref(initalConditions), 
				std::cref(problem),
				std::ref(currentStateVector))));
		}

		//Initalize our status checkup
		bool allowedCheckup = true;

		//Keep running to check the status of everyone
		while (allowedCheckup)
		{
			//Counter to ensure we are all finished
			short int counter = 0;

			//Get the status of the method
			for (map<unsigned int, unique_ptr<SolverIF>>::const_iterator methodItr = allowedMethods.cbegin(); methodItr != allowedMethods.cend(); ++methodItr)
			{
				//Get the params for this method to check its current time
				const OdeSolverParams& currentParams = params.find(methodItr->first)->second;

				//Print our the current percentage done to terminal
				std::cout << std::setprecision(4) << std::setw(2) << "{" << methodItr->first << ":\t" << 100 * std::fabs((currentParams.currentTime - beginTime) / (endTime - beginTime)) << "% Done}" << "\t";

				//Update if we should continue sampling
				counter += static_cast<short int>(currentParams.lastRun);
			}

			//End line
			std::cout << std::endl;

			//check out counter
			if (counter == allowedMethods.size())
			{
				//Break loop
				allowedCheckup = false;
			}

			//Delay our sampling time
			std::this_thread::sleep_for(std::chrono::seconds(2));
		}

		//Join all the threads to get the results
		for (vector<thread>::iterator threadItr = methodThreads.begin(); threadItr != methodThreads.end(); ++threadItr)
		{
			threadItr->join();
		}
	}
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
/// /// Return the results for a known enum Solver Type for the method. 
/// We return the entire vector of solutions.
/// </summary>
/// <param name="methodType"></param>
/// <returns></returns>
const vector<StateVector>& OdeSolver::getResults(SolverIF::SOLVER_TYPES methodType) const
{
	//Get the method Id
	unsigned int methodId = static_cast<unsigned int>(methodType);

	//Check if the method type being asked is in our map
	map<unsigned int, vector<StateVector>>::const_iterator result = resultMap.find(methodId);

	//If result is not found
	if (result == resultMap.cend())
	{
		throw invalid_argument("Invalid Method");
	}

	//Return our vector
	return result->second;
}

/// <summary>
/// Return the results for a known enum value for the method. 
/// We return the entire vector of solutions.
/// </summary>
/// <param name="methodId"></param>
/// <returns></returns>
const vector<StateVector>& OdeSolver::getResults(const unsigned int methodId) const
{
	//Check if the method type being asked is in our map
	map<unsigned int, vector<StateVector>>::const_iterator result = resultMap.find(methodId);

	//If result is not found
	if (result == resultMap.cend())
	{
		throw invalid_argument("Invalid Method");
	}

	//Return our vector
	return result->second;
}

/// <summary>
/// This method will find the "best" / lowest error method and return the vector of all results.
/// </summary>
/// <returns></returns>
const vector<StateVector>& OdeSolver::getResults() const
{
	//check if we even ran any solvers
	if (resultMap.empty())
	{
		throw runtime_error("No method's results saved");
	}

	//Find the best result (smallest error)
	map<unsigned int, vector<StateVector>>::const_iterator bestResult = std::min_element(resultMap.cbegin(), resultMap.cend(),
		[](const std::pair<unsigned int, vector<StateVector>>& leftMap, const std::pair<unsigned int, vector<StateVector>>& rightMap)
		{
			//Get the max total error at the end
			return leftMap.second.back().getParams().totalError < rightMap.second.back().getParams().totalError;
		});
	
	//Check to see if our results are valid
	if (bestResult == resultMap.cend())
	{
		throw runtime_error("Method not found");
	}

	//Return our vector
	return bestResult->second;
}

/// <summary>
/// /// We take in the enum for the method solver type and a desired time we want to find the solution. 
/// We clamp the results if out of bounds.
/// If we want a desired time at a paticular time unsolved, we interpolate the results and build a new state vector.
/// </summary>
/// <param name="methodType"></param>
/// <param name="time"></param>
/// <returns></returns>
const StateVector OdeSolver::getStateAndTime(SolverIF::SOLVER_TYPES methodType, const double time) const 
{
	//Get the method id
	unsigned int methodId = static_cast<unsigned int>(methodType);

	//Call other method to get results for this current method
	return getStateAndTime(methodId, time);

}

/// <summary>
/// We take in the enum value for the method id and a desired time we want to find the solution. 
/// We clamp the results if out of bounds.
/// If we want a desired time at a paticular time unsolved, we interpolate the results and build a new state vector.
/// </summary>
/// <param name="methodId"></param>
/// <param name="time"></param>
/// <returns></returns>
const StateVector OdeSolver::getStateAndTime(const unsigned int methodId, const double time) const 
{
	//Get an iterator to the method
	map<unsigned int, vector<StateVector>>::const_iterator currentMethod = resultMap.find(methodId);

	//Check to see if the method exsits
	if (currentMethod == resultMap.cend())
	{
		throw invalid_argument("Method Invalid");
	}
	//Start searching for the conditions
	else
	{
		//Get a handle to our current results
		const vector<StateVector>& currentResults = currentMethod->second;

		//Check if we need to clamp our results if a time is outside our bounds
		if (currentResults.back().getParams().currentTime <= time)
		{
			return StateVector(currentResults.back().getState(), currentResults.back().getParams());
		}
		else if (currentResults.at(0).getParams().currentTime >= time)
		{
			return StateVector(currentResults.at(0).getState(), currentResults.at(0).getParams());
		}
		//We do not need to clamp
		else
		{
			//This is the iterator to hold the state vector before we pass
			vector<StateVector>::const_iterator beforePassResult;

			//This is the iterator to hold the state vector just after we pass
			vector<StateVector>::const_reverse_iterator afterPassResult;

			//Find the place after after the reqested time
			beforePassResult = std::find_if(currentResults.cbegin(), currentResults.cend(), [&time](const StateVector& s) {return s.getParams().currentTime > time; });

			//Find the place before after the reqested time
			afterPassResult = std::find_if(currentResults.crbegin(), currentResults.crend(), [&time](const StateVector& s) {return s.getParams().currentTime <= time; });

			//Check if we can find our current result
			if (afterPassResult == currentResults.crend() || beforePassResult == currentResults.cend())
			{
				throw invalid_argument("Invalid Time given");
			}
			else
			{
				//Get each found iterators corresponding times
				const double& leftTime = afterPassResult->getParams().currentTime;
				const double& rightTime = beforePassResult->getParams().currentTime;

				//Get each found iterators corresponding state
				const valarray<double>& leftState = afterPassResult->getState();
				const valarray<double>& rightState = beforePassResult->getState();

				//Get the error found on the left and right
				const double& leftError = afterPassResult->getParams().totalError;
				const double& rightError = beforePassResult->getParams().totalError;

				//Get the local truncation error
				const double& leftTrunkError = afterPassResult->getParams().currentError;
				const double& rightTrunkError = beforePassResult->getParams().currentError;

				//Get the runtime
				const double& leftRuntime = afterPassResult->getParams().currentRunTime;
				const double& rightRuntime = beforePassResult->getParams().currentRunTime;

				//Copy over the parameters found
				OdeSolverParams tempParams = afterPassResult->getParams();

				//Copy over the desired time to our parameters
				tempParams.currentTime = time;

				//Build the interpolated state
				valarray<double> intpState = leftState * (1.0 - ((time - leftTime) / (rightTime - leftTime))) +
					rightState * ((time - leftTime) / (rightTime - leftTime));

				//Build the interpolated total error
				tempParams.totalError = leftError * (1.0 - ((time - leftTime) / (rightTime - leftTime))) +
					rightError * ((time - leftTime) / (rightTime - leftTime));

				//Build the interpolated truncation error
				tempParams.currentError = leftTrunkError * (1.0 - ((time - leftTime) / (rightTime - leftTime))) +
					rightTrunkError * ((time - leftTime) / (rightTime - leftTime));

				//Build the interpolated run time
				tempParams.currentRunTime = leftRuntime * (1.0 - ((time - leftTime) / (rightTime - leftTime))) +
					rightRuntime * ((time - leftTime) / (rightTime - leftTime));

				//Update the params and our new state 
				StateVector newState(std::move(intpState), std::move(tempParams));

				//Return our result
				return newState;
			}
		}
	}
}

/// <summary>
/// Find the best solution and return that result
/// </summary>
/// <param name=""></param>
/// <returns></returns>
const StateVector OdeSolver::getStateAndTime(const double time) const
{
	//check if we even ran any solvers
	if (resultMap.empty())
	{
		throw runtime_error("No method's results saved");
	}

	//Find the best result (smallest error)
	map<unsigned int, vector<StateVector>>::const_iterator bestResult = std::min_element(resultMap.cbegin(), resultMap.cend(), 
		[](const std::pair<unsigned int, vector<StateVector>>& leftMap, const std::pair<unsigned int, vector<StateVector>>& rightMap)
		{
			//Get the max total error at the end
			return leftMap.second.back().getParams().totalError < rightMap.second.back().getParams().totalError;
		});

	//Check to see if our results are valid
	if (bestResult == resultMap.cend())
	{
		throw runtime_error("Method not found");
	}

	//Return the result for "best" solver
	return getStateAndTime(bestResult->first, time);
}

/// <summary>
/// Build up all our maps for our methods.
/// If the params given make no sense we throw
/// If no methods were built we throw an invalid argument error.
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
		exit(1);
	}

	//Get the methods allowed
	const methodMap& allowedMethods = methods.getMethodMap();

	//Check if no methods were created
	if (allowedMethods.empty())
	{
		throw runtime_error("No valid methods");
		exit(1);
	}
	else
	{
		//Iterate through the allowed methods to build the parameter tables and results table
		for (methodMap::const_iterator methodItr = allowedMethods.cbegin(); methodItr != allowedMethods.cend(); ++methodItr)
		{
			//This is the ID/enum for the method
			const unsigned int methodId = methodItr->first;

			//Add the method id and parameters to our param map
			params.emplace(methodId, generalParams);

			//set up our result map
			resultMap.emplace(methodId, vector<StateVector>());
		}
	}
}

/// <summary>
/// We set up the time stepping scheme in order to build solutions to the desired time.
/// It calls on methods to run each method and build up our tables.
/// After each time step, we append the solution to our results vector.
/// </summary>
/// <param name="methodId"></param>
/// <param name="currentMethod"></param>
/// <param name="currentParameters"></param>
/// <param name="currentTables"></param>
/// <param name="beginTime"></param>
/// <param name="endTime"></param>
/// <param name="initalConditions"></param>
/// <param name="problem"></param>
/// <param name="results"></param>
void OdeSolver::updateNextTimeStep(const unsigned int methodId, unique_ptr<SolverIF>& currentMethod, OdeSolverParams& currentParameters, 
	Richardson& currentTables, const double beginTime, const double endTime, const valarray<double>& initalConditions, const OdeFunIF* problem, vector<StateVector>& results)
{
	//Get our lock
	mutex lock;

	//Lock the thread
	lock.lock();

	//Get this current methods dt
	double& dt = currentParameters.dt;

	//Unlock
	lock.unlock();

	//Save our current time
	double currentTime = beginTime;

	//Save our currentState
	valarray<double> currentState = initalConditions;

	//Add the current time to our parameters
	currentParameters.currentTime = currentTime;

	//Lock
	lock.lock();

	//Try adding the first state vector to our results
	try
	{
		//Add the first result into results
		results.push_back(std::move(StateVector(currentState, currentParameters)));
	}
	catch (exception& e)
	{
		//Unlock our thread
		lock.unlock();

		//Result error code and exit further processing
		std::cerr << e.what();
		exit(1);
	}

	//Unlock
	lock.unlock();

	while (currentTime < endTime)
	{
		try
		{
			//Lock
			lock.lock();

			//Solver for the next time step for the current method
			currentState = buildSolution(currentMethod, methodId, currentTables, currentParameters, currentState, problem, currentTime, endTime);

			//Unlock
			lock.unlock();
		}
		catch (exception& e)
		{
			//Unlock
			lock.unlock();

			//Print our error and exit the program
			cerr << e.what();
			exit(1);
		}

		//Lock
		lock.lock();

		//Update the time
		currentTime += dt;

		//Add the new time to our parameters
		currentParameters.currentTime = currentTime;

		//Try to append the current results to our results map
		try
		{
			//Push back the result
			results.push_back(std::move(StateVector(currentState, currentParameters)));
		}
		catch (exception& e)
		{
			//Unlock
			lock.unlock();

			cerr << e.what();
			exit(1);
		}

		//Unlock
		lock.unlock();
	}
}

const bool OdeSolver::isExplict(const unsigned int methodId) const
{
	switch (methodId)
	{
	case 10:
	case 20:
	case 30:
	{
		return true;
	}
	case 40:
	case 50:
	{
		return false;
	}
	default:
	{
		throw std::invalid_argument("Method Id not found");
	}
	}
}