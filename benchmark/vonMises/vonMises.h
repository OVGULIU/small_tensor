#ifndef OPENUS_FE_MATERIAL_vonMises_H_
#define OPENUS_FE_MATERIAL_vonMises_H_

#include <smalltensor/smalltensor.h>

namespace openus{

class vonMises
{
typedef smalltensor::tensor2<float,3,3> Mat33 ; 
typedef smalltensor::tensor4<float,3,3,3,3> Tensor4 ;

public:
	vonMises();
	void Initialize() ; 
	void Update() ; 

	Tensor4 const& GetStiffnessTensor() const ;

	Tensor4 _stiff_tensor ;

	Mat33 _iter_stress ;
	Mat33 _iter_strain ;
	Mat33 _iter_plastic_strain ;
	Mat33 _iter_back_stress ;

	Mat33 _commit_stress ;
	Mat33 _commit_strain ;
	Mat33 _commit_plastic_strain ;
	Mat33 _commit_back_stress ;

	static const Mat33 _zero_strain ;
	static const Mat33 _zero_stress ;
	static const Mat33 kronecker_delta ;

	smalltensor::eindex<'I'> I;
	smalltensor::eindex<'J'> J;
	smalltensor::eindex<'K'> K;
	smalltensor::eindex<'L'> L;
	smalltensor::eindex<'P'> P;
	smalltensor::eindex<'Q'> Q;

};


}


#endif

