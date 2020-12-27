#include "Richardson.h"

void Richardson::BuildTables(const size_t tableSize)
{
	//Save the table size
	N = static_cast<unsigned int>(tableSize);

	//Resize the table
	result.resize(tableSize);

	//Resize each column
	for (int i = 0; i < result.size(); ++i)
	{
		result[i].resize(tableSize);
	}
}

void Richardson::initalizeSteps(const double& reduct, const double& dt)
{
	//Save the reduction factor
	reductionFactor = reduct;

	//Save dt
	stepSize = dt;
}

void Richardson::operator()(const size_t rowIndx, const size_t colIndx, crvec currentResult)
{
	try
	{
		result[rowIndx][colIndx] = currentResult;
	}
	catch (exception& e)
	{
		cerr << e.what();
	}
}