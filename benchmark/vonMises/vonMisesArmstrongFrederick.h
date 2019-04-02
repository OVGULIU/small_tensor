#ifndef OPENUS_FE_MATERIAL_vonMisesArmstrongFrederick_H_
#define OPENUS_FE_MATERIAL_vonMisesArmstrongFrederick_H_

#include <openus/fe/material/FEMaterialBase.h>

namespace openus{


class vonMisesArmstrongFrederick: public FEMaterialBase
{
public:
	vonMisesArmstrongFrederick(entt::registry<>* registry);
	void Initialize() ; 
	void Update() ; 

	math::Tensor4f const& GetStiffnessTensor() const ;

	math::Tensor4f _stiff_tensor ;

};


}


#endif

