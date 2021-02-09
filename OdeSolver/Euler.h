#pragma once

#include <valarray>

#include "OdeFunIF.h"
#include "SolverIF.h"

//Convience for writing out methods
using std::valarray;
using vec = valarray<double>;
using crvec = const vec&;
using rvec = vec&;

// Class derrived from the SolverIF to support the Euler time stepping scheme
class Euler : public SolverIF
{
protected:

	// Hold the functions derivative vector at the current time step. Used to "move" to current state vector in time
	vec k1;

private:

	// Update the vectors that are used to appoximate the function vectors derivative at other time steps
	virtual void initalizeSolverVectors() override;

public:

	// Using default constructor
	Euler() = default;

	// Using default copy constructor
	Euler(const Euler&) = default;

	// Using default assignment operator
	Euler& operator=(const Euler&) = default;

	// Using default destructor
	virtual ~Euler() = default;

	// Initalize the current state vector and solving helper vectors
	virtual void initalize(crvec) override;

	/// Update the current vector's state for explct methods
	virtual rvec update(crvec, rvec, const double&, const double&, const int&, const OdeFunIF*) override;

	//Get the next time step for rvec for implict methods
	virtual rvec update(crvec, rvec, const double&, const double&, const int&, const OdeFunIF*, const double&, const double&) override;

	/// Return the error order of the Euler method
	virtual const double getErrorOrder() const override;
};

