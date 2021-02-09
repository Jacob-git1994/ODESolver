#pragma once

#include <valarray>

#include "Euler.h"
#include "OdeFunIF.h"

using std::valarray;

class RK2 : public Euler
{
protected:

	//Vector to hold function vector
	valarray<double> k2;

private:

	//Initalize our solving vectors
	virtual void initalizeSolverVectors() override;

public:

	//Default constructor
	RK2() = default;

	//Default copy constructor
	RK2(const RK2&) = default;

	//Default Assign operator
	RK2& operator=(const RK2&) = default;

	//Default destructor
	virtual ~RK2() = default;

	//Initalize the vector
	virtual void initalize(const valarray<double>&) override;

	//Get the next time step for rvec for explict methods
	virtual rvec update(const valarray<double>&, valarray<double>&, const double&, const double&, const int&, const OdeFunIF*) override;

	//Get the next time step for rvec for implict methods
	virtual rvec update(crvec, rvec, const double&, const double&, const int&, const OdeFunIF*, const double&, const double&) override;

	//Get the power of the error
	virtual const double getErrorOrder() const override;
};

