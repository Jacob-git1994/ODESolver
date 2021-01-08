#pragma once

#include <array>

#include "SolverIF.h"

using std::array;

class OdeSolverParams
{
private:

	//Check user inputs
	bool checkUSerInputs() const;

public:

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

	//Richardson table generation
	size_t minTableSize;
	size_t maxTableSize;
	size_t currentTableSize;

	//Problem Specifics
	bool isStiff; //Stiff PDEs
	bool isLarge; //The problem requires a large table
	bool isFast; //The problem evolves quickly

	//Construtors
	OdeSolverParams(const array<bool,5>, const array<double,2>, const array<size_t,2>, const array<size_t,3>);

	//Copy Constructor
	OdeSolverParams(const OdeSolverParams&);

	//Assign Operator
	const OdeSolverParams& operator=(const OdeSolverParams&);
};

OdeSolverParams::OdeSolverParams(const array<bool, 5> allowedMethods = {true,false,false,false,false},
	const array<double, 2> errorBounds = {.001,.01},
	const array<size_t, 2> richLevelBounds = {4,8},
	const array<size_t, 3> problemSpecifics = {false,false,false}) :
	useEuler(allowedMethods[0]),
	useRK2(allowedMethods[1]),
	useRK4(allowedMethods[2]),
	useImplictEuler(allowedMethods[3]),
	useCrank(allowedMethods[4]),
	lowerError(errorBounds[0]),
	upperError(errorBounds[1]),
	minTableSize(richLevelBounds[0]),
	maxTableSize(richLevelBounds[1]),
	isStiff(problemSpecifics[0]),
	isLarge(problemSpecifics[1]),
	isFast(problemSpecifics[2]),
	currentError(0.0),
	currentTableSize(1)
{
	//Do nothing here
}

