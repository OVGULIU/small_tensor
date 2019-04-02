#ifndef OPENUS_FE_MATERIAL_vonMisesLinearHardening_H_
#define OPENUS_FE_MATERIAL_vonMisesLinearHardening_H_

// #include <openus/fe/material/FEMaterialBase.h>
#include "vonMises.h"

namespace openus{

class vonMisesLinearHardening: public vonMises
{
public:
	vonMisesLinearHardening();
	// void Initialize() ; 
	// void Update() ; 

	Mat33 kinematic_derivative(Mat33 const& m, Mat33 const& back_stress) const override final;

	float _kinematic_harden_rate = 4 * 10 ;

};


}


#endif

