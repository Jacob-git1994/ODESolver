#include "Richardson.h"

void Richardson::BuildTables(const size_t tableSize, const size_t vecSize)
{
	try 
	{
		//Resize the table
		result.resize(tableSize);

		//Resize each column
		for (size_t i = 0; i < result.size(); ++i)
		{
			result[i].resize(tableSize);
		}

		//Initalize all the stored vectors
		for (size_t i = 0; i < result.size(); ++i)
		{
			for (size_t j = 0; j < result.size(); ++j)
			{
				result[i][j].resize(vecSize);
			}
		}
	}
	catch (exception& e)
	{
		cerr << e.what();
		exit(-1);
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
		exit(-1);
	}
}

double Richardson::normedError() const
{
	//Get the error vector
	vec error = result[result.size() - 1][result.size() - 1] - result[result.size() - 2][result.size() - 2];

	//Calculate the norm
	double normVal = 0.;
	for (size_t i = 0; i < error.size(); ++i)
	{
		normVal += error[i] * error[i];
	}

	//Return the sqrt of the norm
	return sqrt(normVal);
}

double Richardson::error(rvec bestResult)
{
	//Get ref to the current object
	Richardson& currentTable = *this;

	//Iterate through the rows of the table
	for (size_t i = 0; i < result.size()-1; ++i)
	{
		//Iterate through the columns
		for (size_t j = 0; j <= i; ++j)
		{
			//Get the updated result
			vec updatedResult = (pow(reductionFactor, static_cast<double>(j) + 1.)
				* result[i+1][j] - result[i][j])
				/ (pow(reductionFactor, static_cast<double>(j) + 1.) - 1.);

			//Save the updated result to the table
			currentTable(i + 1, j + 1, updatedResult);
		}
	}

	//Get the last result as that is the "best one"
	bestResult = result[result.size() - 1][result.size() - 1];
	currentNormError = normedError();

	//return error
	return currentNormError;
}