#pragma once
#include "LinAlgHelperBase.h"
class FirstOrderScheme : public LinAlgHelperBase
{
private:

	//Override our function for the jacobian
	virtual void getJacobian(const OdeFunIF*, const double&, const double&) override;
};

