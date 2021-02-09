#pragma once

#include <valarray>

#include "Euler.h"
#include "OdeFunIF.h"

using std::valarray;

class RK4 : public Euler
{
protected:

	//Vector to hold the function vector
	vec k2;
	vec k3;
	vec k4;

private:
	//Initalize solving vectors
	virtual void initalizeSolverVectors() override;

public:

	//Default Constrcutor
	RK4() = default;

	//Default Copy Constructor
	RK4(const RK4&) = default;

	//Default Assign operator
	RK4& operator=(const RK4&) = default;

	//Default Destrutor
	virtual ~RK4() = default;

	//Initalize the vector
	virtual void initalize(crvec) override;

	//Get the next time step for rvec for explict methods
	virtual rvec update(crvec, rvec, const double&, const double&, const int&, const OdeFunIF*) override;

	//Get the next time step for rvec for implict methods
	virtual rvec update(crvec, rvec, const double&, const double&, const int&, const OdeFunIF*, const double&, const double&) override;

	//Get the power of the error
	virtual const double getErrorOrder() const override;
};

