#pragma once

#include "Euler.h"

#include "MethodWrapperBase.h"
class MethodEuler : virtual public MethodWrapperBase
{
protected:

	//Build our solvers
	virtual void buildSolvers() override;

public:

	//Default Constrcutor
	MethodEuler() = default;

	//Delete the copy constructor
	MethodEuler(const MethodEuler&) = delete;

	//Default Deconstrcutor
	virtual ~MethodEuler() = default;
};

