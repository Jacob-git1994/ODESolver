#pragma once
class OdeSolverParams
{
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
};

