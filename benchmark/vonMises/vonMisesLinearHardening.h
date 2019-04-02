#ifndef OPENUS_FE_MATERIAL_vonMisesLinearHardening_H_
#define OPENUS_FE_MATERIAL_vonMisesLinearHardening_H_

#include <openus/fe/material/FEMaterialBase.h>

namespace openus{


class vonMisesLinearHardening: public FEMaterialBase
{
public:
	vonMisesLinearHardening(entt::registry<>* registry);
	void Initialize() ; 
	void Update() ; 

	math::Tensor4f const& GetStiffnessTensor() const ;

	math::Tensor4f _stiff_tensor ;

};


}


#endif

