#ifndef OPENUS_FE_MATERIAL_vonMisesArmstrongFrederick_H_
#define OPENUS_FE_MATERIAL_vonMisesArmstrongFrederick_H_

// #include <openus/fe/material/FEMaterialBase.h>
#include "vonMises.h"

namespace openus{

class vonMisesArmstrongFrederick: public vonMises
{
public:
	vonMisesArmstrongFrederick();
	// void Initialize() ; 
	// void Update() ; 

	Mat33 kinematic_derivative(Mat33 const& m, Mat33 const& back_stress) const override final;

	float _AF_ha = 3e3 ; 
	float _AF_cr = 500 ;	
};


}


#endif

