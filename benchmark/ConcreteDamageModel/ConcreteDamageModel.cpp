
#include "ConcreteDamageModel.h"

using namespace smalltensor ; 

const tensor2<float,3,3> ConcreteDamageModel::kronecker_delta("identity");


ConcreteDamageModel::ConcreteDamageModel(
		float E_in,
		float v_in,
		float beta_in,
		float K_in, 
		float r_0_tensile_in,
		float r_0_compress_in,
		float mesh_A_tensile_in,
		float mesh_A_compress_in,
		float mesh_B_compress_in
	)
{
	_E = E_in ; 
	_v = v_in ; 
	_beta_plastic_intensity = beta_in ; 
	_K_compress = K_in ; 
	_r_0_tensile = r_0_tensile_in ; 
	_r_0_compress = r_0_compress_in ; 

	_r_tensile = _r_0_tensile ; 
	_r_compress = _r_0_compress ; 

	_mesh_A_tensile = mesh_A_tensile_in ; 
	_mesh_A_compress = mesh_A_compress_in ; 
	_mesh_B_compress = mesh_B_compress_in ; 

	_d_tensile = 0.f ; 
	_d_compress = 0.f ; 


	double lambda = ( _v * _E ) / ( ( 1 + _v ) * ( 1 - 2 * _v ) );
	double mu = _E / ( 2 * ( 1 + _v ) );
	_Ee(I,J,K,L) = mu * (
					  kronecker_delta(I,K) * kronecker_delta(J,L) 
					+ kronecker_delta(I,L) * kronecker_delta(J,K)
					)
	                + lambda * kronecker_delta(K,L)*kronecker_delta(I,J) ;
    _D_inv(I,J,K,L) = 1 / _E * (
    		(1 + _v)/2 * 
    			( kronecker_delta(I,K) * kronecker_delta(J,L) 
				+ kronecker_delta(I,L) * kronecker_delta(J,K)
				)
			- _v * kronecker_delta(I,J) * kronecker_delta(K,L) 
    	) ; 
}
// ~ConcreteDamageModel();


int 
ConcreteDamageModel::setTrialStrainIncr(tensor2<float,3,3> const& strain_incr){
	_trial_strain(I,J) = _trial_strain(I,J) + strain_incr(I,J) ; 

	// Algorithm Box-1 Computer Effective Stress
	computeEffectiveStress(strain_incr);

	// Algorithm Box-2 Computer Cauchy Stress
	computeCauchyStress();

}


// Algorithm Box-1 Computer Effective Stress
int 
ConcreteDamageModel::computeEffectiveStress(tensor2<float,3,3> const& strain_incr){
	// Step (i): compute the trial stress.
	_effective_stress_trial(I,J) = _effective_stress_commit(I,J)
	 		+ _Ee(I,J,K,L) * strain_incr(K,L) ; 

	// Step (ii): is beta zero?
	if ( std::fabs(_beta_plastic_intensity) < 1e-10 ){
		// Step (ii):  No plasticity.
		_effective_stress_next = _effective_stress_trial ;
		return 1 ; 
	}

	// Step (iii): split stress to tensile and compressive
	SplitStress(_effective_stress_trial, _effective_stress_tensile, _effective_stress_compress)

	const float compress_equiv_stress = GetEquivStressCompress(_effective_stress_compress, _K_compress) ; 

	if ( compress_equiv_stress < _r_compress ){
		// Step (iii):  No evolution for _d_compress and plastic_strain.
		_effective_stress_next = _effective_stress_trial ;
		return 1; 
	}
	
	// Step (iv): plastic evolution? 
	const float trial_stress_norm = std::sqrt( _effective_stress_trial(I,J) * _effective_stress_trial(I,J) ); 

	tensor2<float,3,3> unit_trial_stress ; 
	unit_trial_stress(I,J) = _effective_stress_trial(I,J) / trial_stress_norm ; 

    const float stress_strain_contract = unit_trial_stress(I,J) * strain_incr(I,J) ; 
    if ( stress_strain_contract <= 0 ){
    	// Step (iv): No plastic evolution.
    	_effective_stress_next = _effective_stress_trial ;
    	return 1; 
    }
    const float lambda = 1 - _beta_plastic_intensity / trial_stress_norm * _E 
                  * Heaviside( _d_compress /*incremental?*/ ) * Macaulay(stress_strain_contract)  ; 
    tensor2<float,3,3> effective_stress_hat ;
    tensor2<float,3,3> effective_stress_hat_tensile ;
    tensor2<float,3,3> effective_stress_hat_compress ;
    effective_stress_hat(I,J) = lambda * _effective_stress_trial(I,J) ; 
    SplitStress(effective_stress_hat, effective_stress_hat_tensile, effective_stress_hat_compress) ; 
    const float compress_equiv_stress_hat = GetEquivStressCompress(effective_stress_hat_compress) ;
    if ( compress_equiv_stress_hat < _r_compress){
    	// Step (iv): No evolution for _d_compress and plastic_strain.
    	_effective_stress_next = _effective_stress_trial ;
    	return 1 ;
    }

    // TODO: Need evolution of plastic_strain.
    _effective_stress_next = effective_stress_hat ; 

}

int 
ConcreteDamageModel::computeCauchyStress(){
	tensor2<float,3,3> effective_stress_tensile ; 
	tensor2<float,3,3> effective_stress_compress ; 
	SplitStress(_effective_stress_next, effective_stress_tensile , effective_stress_compress) ; 
	const float equiv_stress_compress = GetEquivStressCompress(effective_stress_compress, _K_compress) ; 
	const float equiv_stress_tensile  = std::sqrt( effective_stress_tensile(I,J) * 
					_D_inv(I,J,K,L) * effective_stress_tensile(K,L) ) ; 

	_r_compress = std::max(_r_compress, equiv_stress_compress) ; 
	_r_tensile = std::max(_r_tensile, equiv_stress_tensile) ; 

	updateDamageVariables();

	_cauchy_stress(I,J) = ( 1 - _d_tensile ) * effective_stress_tensile(I,J)
						+ ( 1 + _d_compress) * effective_stress_compress(I,J) ; 
}

int ConcreteDamageModel::updateDamageVariables(){
	// TODO: Need the iterative of damage_variables for iterative FEM.

	if (_r_tensile >= _r_0_tensile) {
		_d_tensile = 1 - _r_0_tensile / _r_tensile * 
			std::exp(_mesh_A_tensile * (1 - _r_tensile / _r_0_tensile)) ; 
	}

	if (_r_compress >= _r_0_compress) {
		_d_compress = 1 - _r_0_compress / _r_compress * ( 1 - _mesh_A_compress)
		  	- _mesh_A_compress * std::exp(_mesh_B_compress * (1 - _r_compress / _r_0_compress) )  ;
	}
}

void 
ConcreteDamageModel::commitState(){

	_commit_strain = _trial_strain ; 

}

smalltensor::tensor2<float,3,3> 
ConcreteDamageModel::getStressTensor() const{
	return _cauchy_stress ; 
}

smalltensor::tensor2<float,3,3> 
ConcreteDamageModel::getStrainTensor() const{
	return _trial_strain ; 
}


