#pragma once

#include <array>
#include <stdexcept>

#include "SolverIF.h"

using std::array;
using std::invalid_argument;

class OdeSolverParams
{

public:

	//Check user inputs
	inline bool checkUserInputs() const;

	//Allowed Methods
	bool useEuler;
	bool useRK2;
	bool useRK4;
	bool useImplictEuler;
	bool useCrank;

	//Error Bounds Allowed
	double lowerError;
	double upperError;
	double currentError;
	bool satifiesError;

	//Allowed deleta time
	double minDt;
	double maxDt;
	double dt;

	//Error of the constant
	double c;

	//Parameters for dt
	double isDtClamped;

	//Time constraints
	double currentRunTime;

	//Richardson table generation
	size_t minTableSize;
	size_t maxTableSize;
	size_t currentTableSize;

	//Richardson reduction factors
	double redutionFactor;

	//Problem Specifics
	bool isStiff; //Stiff PDEs
	bool isLarge; //The problem requires a large table
	bool isFast; //The problem evolves quickly

	//Construtors
	inline OdeSolverParams(const array<bool, 5>&, const array<double, 2>&, const array<double, 2>&, const array<size_t, 2>&, const array<size_t, 3>&, const double&);

	//Copy Constructor
	inline OdeSolverParams(const OdeSolverParams&) = default;

	//Assign Operator
	inline const OdeSolverParams& operator=(const OdeSolverParams&);
};

OdeSolverParams::OdeSolverParams(const array<bool, 5>& allowedMethods = {true,false,false,false,false},
	const array<double, 2>& errorBounds = {.0001,.001},
	const array<double, 2>& dtBounds = {.01,.1},
	const array<size_t, 2>& richLevelBounds = {4,8},
	const array<size_t, 3>& problemSpecifics = {false,false,false},
	const double& reductionFactorIn = 2.) :
	useEuler(allowedMethods[0]),
	useRK2(allowedMethods[1]),
	useRK4(allowedMethods[2]),
	useImplictEuler(allowedMethods[3]),
	useCrank(allowedMethods[4]),
	lowerError(errorBounds[0]),
	upperError(errorBounds[1]),
	minDt(dtBounds[0]),
	maxDt(dtBounds[1]),
	minTableSize(richLevelBounds[0]),
	maxTableSize(richLevelBounds[1]),
	isStiff(problemSpecifics[0]),
	isLarge(problemSpecifics[1]),
	isFast(problemSpecifics[2]),
	currentError(0.0),
	currentTableSize(1),
	currentRunTime(0.0),
	dt(.01),
	redutionFactor(reductionFactorIn),
	isDtClamped(false),
	satifiesError(true),
	c(-1.0)
{
	//If the inputs are invalid we do no want to continue
	if (!checkUserInputs())
	{
		throw invalid_argument("Invalid ODE parameters");
	}
}

bool OdeSolverParams::checkUserInputs() const
{
	//Initalize our arguments
	bool goodArgs = true;

	//Make sure the error inputs are valid
	if (!((lowerError > 0. && isfinite(lowerError) && lowerError < upperError) && (upperError > 0 && isfinite(upperError) && upperError > lowerError)))
	{
		goodArgs = false;
	}

	//Make sure the table arguments are valid
	if (!(minTableSize > 1 and maxTableSize > 2 && minTableSize < maxTableSize))
	{
		goodArgs = false;
	}

	//Make sure the run times are valid
	if (!((minDt > 0. && isfinite(minDt) && minDt < maxDt) && (maxDt > 0 && isfinite(maxDt) && maxDt > minDt)))
	{
		goodArgs = false;
	}

	if (!(redutionFactor > 1))
	{
		goodArgs = false;
	}

	//Return if the arguments are valid
	return goodArgs;
}

const OdeSolverParams& OdeSolverParams::operator=(const OdeSolverParams& params)
{
	//Copy all the parameters over
	useEuler = params.useEuler;
	useRK2 = params.useRK2;
	useRK4 = params.useRK4;
	useImplictEuler = params.useImplictEuler;
	useCrank = params.useCrank;
	lowerError = params.lowerError;
	upperError = params.upperError;
	currentError = params.currentError;
	minDt = params.minDt;
	maxDt = params.maxDt;
	currentRunTime = params.currentRunTime;
	minTableSize = params.minTableSize;
	maxTableSize = params.maxTableSize;
	currentTableSize = params.currentTableSize;
	isStiff = params.isStiff; 
	isLarge = params.isLarge;
	isFast = params.isFast; 
	dt = params.dt;
	redutionFactor = params.redutionFactor;
	isDtClamped = params.isDtClamped;
	satifiesError = params.satifiesError;
	c = params.c;

	//Return this
	return *this;
}