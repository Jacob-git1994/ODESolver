#pragma once

#include "MethodEuler.h"
#include "MethodRK4.h"
#include "MethodWrapperBase.h"

class MethodEulerRK4 : public MethodEuler, public MethodRK4
{
protected:

	//Build our solver
	virtual void buildSolvers() override;

public:

	//Default Constrcutor
	MethodEulerRK4() = default;

	//Delete the copy constructor
	MethodEulerRK4(const MethodEulerRK4&) = delete;

	//Default Deconstrcutor
	virtual ~MethodEulerRK4() = default;
};

