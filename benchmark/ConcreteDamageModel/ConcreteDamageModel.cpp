


#include "ConcreteDamageModel.h"

using namespace smalltensor ; 

const tensor2<float,3,3> ConcreteDamageModel::kronecker_delta("identity");


ConcreteDamageModel::ConcreteDamageModel(

	)
{
	_d_tensile = 0.f ; 
	_d_compress = 0.f ; 
	_r_0_tensile = 0.f ; 
	_r_0_compress = 0.f ; 
	_r_tensile = 0.f ; 
	_r_compress = 0.f ; 

	double lambda = ( _v * _E ) / ( ( 1 + _v ) * ( 1 - 2 * _v ) );
	double mu = _E / ( 2 * ( 1 + _v ) );
	_Ee(I,J,K,L) = mu * (
					  kronecker_delta(I,K) * kronecker_delta(J,L) 
					+ kronecker_delta(I,L) * kronecker_delta(J,K)
					)
	                + lambda * kronecker_delta(K,L)*kronecker_delta(I,J) ;
}
// ~ConcreteDamageModel();


int 
ConcreteDamageModel::setTrialStrainIncr(tensor2<float,3,3> const& strain_incr){
	_effective_stress_trial(I,J) = _effective_stress_commit(I,J)
	 		+ _Ee(I,J,K,L) * strain_incr(K,L) ; 

	compute_effective_stress();

	compute_cauchy_stress();

}

int 
ConcreteDamageModel::compute_effective_stress(){
	split_
}

int 
ConcreteDamageModel::compute_cauchy_stress(){

}

void 
ConcreteDamageModel::commitState(){

}

smalltensor::tensor2<float,3,3> 
ConcreteDamageModel::getStressTensor() const{

}

smalltensor::tensor2<float,3,3> 
ConcreteDamageModel::getStrainTensor() const{

}


