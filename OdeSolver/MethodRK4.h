#pragma once

#include "MethodWrapperBase.h"
#include "RK4.h"

class MethodRK4 : virtual public MethodWrapperBase
{
protected:

	//Build our solvers
	virtual void buildSolvers() override;

public:

	//Default Constrcutor
	MethodRK4() = default;

	//Delete the copy constructor
	MethodRK4(const MethodRK4&) = delete;

	//Default Deconstrcutor
	virtual ~MethodRK4() = default;
};

