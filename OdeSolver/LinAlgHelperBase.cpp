#include "LinAlgHelperBase.h"

void LinAlgHelperBase::getFuncDer(const OdeFunIF* problemIn, const double& currentTime)
{
	//Update our function vector with the current state
	problemIn->operator()(funcVec, guessLeft, currentTime);
}

const valarray<double> LinAlgHelperBase::solveSystem()
{
	//Our result vector
	valarray<double> result = funcVec;

	//Matrix size (assuming it is square..will need a check here)
	const size_t jSize = A.size();

	//Simple inefficent soler until later ... just need it working .. lol
	for (size_t i = 0; i < jSize; ++i)
	{
		//Check if 0
		if (A[i][i] == 0.0)
		{
			throw std::runtime_error("Division by Zero");
		}

		for (size_t j = 0; j < jSize; ++j)
		{
			double ratio = A[j][i] / A[i][i];

			for (size_t k = 0; k < jSize; ++k)
			{
				A[j][k] += -ratio * A[i][k];
			}
		}
	}

	//Back subsitution end
	result[jSize - 1] = A[jSize - 2][jSize - 1] / A[jSize - 1][jSize - 1];

	//Back subsitution
	for (size_t i = jSize - 2; i >= 0; --i)
	{
		for (size_t j = i + 1; j < jSize; ++j)
		{
			result[i] += -A[i][j] * result[j];
		}

		result[i] *= (1.0 / A[i][i]);
	}

	//Return our result
	return result;
}

LinAlgHelperBase::LinAlgHelperBase(const double& errorTolIn, const double& dtIn, const unsigned int& maxItrIn) :
	errorTol(errorTolIn),
	dt(dtIn),
	maxIter(maxItrIn)
{
	//Nothing else to do here
}

const valarray<double>& LinAlgHelperBase::solve(const double& currentTime, const double& methodDt, const valarray<double>& currentState, const OdeFunIF* problemIn)
{
	//Generate our first pair of guesses
	guessLeft = currentState;
	guessRight = currentState + methodDt * problemIn->operator()(guessRight, guessLeft, currentTime);

	//Generate our estimated error
	double error = 99999;

	//Set our current iteration count
	unsigned int iter = 0;

	//Initalize our error vector
	valarray<double> errorVec;

	valarray<double> storageVec;

	//Iterate
	do
	{
		//Generate the Jacobian for A
		getJacobian(problemIn, currentTime, methodDt);

		//Generate the function derv vector
		getFuncDer(problemIn, currentTime);

		//Update our function vector
		funcVec = guessLeft - currentState - methodDt * problemIn->operator()(storageVec, guessLeft, currentTime + methodDt);

		try
		{
			//Solve our system
			guessLeft += -solveSystem();
		}
		catch (std::exception& e)
		{
			std::cerr << e.what();
		}

		//Find the error vector
		errorVec = guessLeft - guessRight;

		//Get the absolute value of each element
		std::for_each(std::begin(errorVec), std::end(errorVec), 
			[](double& a)
			{
				a *= a;
			});

		//Get the 2 normed error
		error = std::sqrt(errorVec.sum());

		//Update the vectors
		guessRight = guessLeft;

		//Break once these conditions occur
	} while (error < errorTol && iter < ++maxIter);

	//Return our best iteration
	return guessLeft;
}
