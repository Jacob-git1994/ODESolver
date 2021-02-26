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

	isBuilt = true;
}

void Richardson::initalizeSteps(const double& reduct, const double& dt)
{
	//Save the reduction factor
	reductionFactor = reduct;

	//Save dt
	stepSize = dt;
}

void Richardson::append(const size_t rowIndx, const size_t colIndx, valarray<double>&& currentResult)
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

void Richardson::append(const size_t rowIndx, const size_t colIndx, const valarray<double>& currentResult)
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
	//vec error = result[result.size() - 1][result.size() - 1] - result[0][0];

	for_each(std::begin(error), std::end(error), [](double& elIn)
		{
			elIn = std::abs(elIn);
		});

	return error.max();
}

const double Richardson::error(rvec bestResult, double& c)
{
	//Get ref to the current object
	Richardson& currentTable = *this;

	//Iterate through the rows of the table
	for (size_t i = 1; i < result.size(); ++i)
	{
		//Iterate through the columns
		for (size_t j = 0; j < i; ++j)
		{
			//Get the updated result
			vec updatedResult = (pow(reductionFactor, static_cast<double>(j) + 1.)
				* result[i][j] - result[i - 1][j])
				/ (pow(reductionFactor, static_cast<double>(j) + 1.) - 1.);

			//Save the updated result to the table
			currentTable.append(i, j + 1, std::move(updatedResult));
		}
	}

	//Get the last result as that is the "best one"
	bestResult = result[result.size() - 1][result.size() - 1];
	currentNormError = normedError();

	//Set c to our approximaation of convergence
	c = abs(log(currentNormError) / log(stepSize));

	//return error
	return currentNormError;
}

const size_t Richardson::getTableSize() const
{
	return result.size();
}

const double Richardson::getReductionFactor() const
{
	return reductionFactor;
}
