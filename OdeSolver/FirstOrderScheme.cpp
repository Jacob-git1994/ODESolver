#include "FirstOrderScheme.h"

void FirstOrderScheme::getJacobian(const OdeFunIF* problemIn, const double& currentTime, const double& methodDt)
{
	valarray<double> temp1, temp2;

	//Initalize our Jacobian
	for (size_t i = 0; i < A.size(); ++i)
	{
		//Copy over the guess for modification
		valarray<double> copyGuess = guessLeft;

		//Update the element we are taking the derivative for
		copyGuess[i] += dt;

		A[i] = -(methodDt / dt) * (problemIn->operator()(temp1, copyGuess, currentTime + methodDt) - (problemIn->operator()(temp2, guessLeft, currentTime + methodDt)));
	}

	//Add one down the eyes
	for (size_t i = 0; i < A.size(); ++i)
	{
		A[i][i] += 1.0;
	}

}
