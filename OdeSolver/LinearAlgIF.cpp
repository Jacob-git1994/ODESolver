#include "LinearAlgIF.h"

/// <summary>
/// We will build the derivative function for the current state and time input
/// </summary>
/// <param name=""></param>
/// <param name=""></param>
/// <param name=""></param>
/// <param name=""></param>
void LinearAlgIF::buildFuncDer(const OdeFunIF* problemIn, valarray<double>& resultIn, const valarray<double>& currentStateIn, const double& currentTimeIn)
{
	//Generate the result
	problemIn->operator()(resultIn, currentStateIn, currentTimeIn);
}

/// <summary>
/// Solve the system of linear equations. Using a simple model for now - later parallize it
/// </summary>
/// <returns></returns>
const valarray<double> LinearAlgIF::solveSystem()
{
	//Vector size
	const size_t vectorSize = leftDfDt.size();

	//Temp storage for the resukt
	valarray<double> result(vectorSize);

	//Iterate through the elements and solve // Will need to think harder here
	for (size_t i = 0; i < vectorSize; ++i)
	{
		for (size_t j = 0; j < vectorSize; ++j)
		{
			for (size_t k = 0; k < vectorSize; ++k)
			{

			}
		}
	}
	return result;
}
/// <summary>
/// Initalize our variables and sizes for future uses
/// </summary>
/// <param name="vectorSizeIn"></param>
/// <param name="allowedErrorIn"></param>
void LinearAlgIF::initalize(const size_t& vectorSizeIn, const double& allowedErrorIn)
{
	//Set the allowed error
	allowedError = allowedErrorIn;

	//Initalize the matrix and vectors with the vector size
	leftDfDt.resize(vectorSizeIn);
	rightDfDt.resize(vectorSizeIn);
	rightResult.resize(vectorSizeIn);
	leftResult.resize(vectorSizeIn);
	J.resize(vectorSizeIn);
	for_each(std::begin(J), std::end(J), [&vectorSizeIn](valarray<double>& rowIn)
		{
			rowIn.resize(vectorSizeIn);
		});
	
}
