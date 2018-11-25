
#include "../../smalltensor/smalltensor.h"
#include "ConcreteDamageModel_utils.hpp"
#include <limits>
#include <vector>
#include <cmath>

using namespace smalltensor ; 
class ConcreteDamageModel
{
public:
	ConcreteDamageModel(
		float E_in,
		float v_in,
		float beta_in,
		float K_in, 
		float r_0_tensile_in,
		float r_0_compress_in,
		float mesh_A_tensile_in,
		float mesh_A_compress_in,
		float mesh_B_compress_in
	);
	// ~ConcreteDamageModel();
	

	int setTrialStrainIncr(tensor2<float,3,3> const& strain_incr);
	void commitState();
	tensor2<float,3,3> getStressTensor() const;
	tensor2<float,3,3> getStrainTensor() const;


private:
	int computeEffectiveStress(tensor2<float,3,3> const& strain_incr) ; 
	int computeCauchyStress() ; 
	int updateDamageVariables() ; 
private:
	float _d_tensile ; 
	float _d_compress ; 
	
	float _r_tensile ; 
	float _r_compress ; 

	float _E;                              // Elastic Modulus
	float _v;                              // Poisson Ratio
	float _rho;                            // Density
	float _K_compress;                     // K compressive
	float _r_0_tensile ;                   // Tensile Damage Radius
	float _r_0_compress ;                  // Compress Damage Radius
	float _mesh_A_tensile ;                // Unique Parameter for mesh-objectivity.
	float _mesh_A_compress ;               // Unique Parameter for mesh-objectivity.
	float _mesh_B_compress ;               // Unique Parameter for mesh-objectivity.

	float _beta_plastic_intensity ;

	tensor2<float,3,3> _effective_stress ;                    // 

	tensor2<float,3,3> _effective_stress_trial ;
	tensor2<float,3,3> _effective_stress_next ;


	tensor2<float,3,3> _cauchy_stress ;


	tensor2<float,3,3> _effective_stress_compress;            // -
	tensor2<float,3,3> _effective_stress_tensile;             // +


	tensor2<float,3,3> _trial_strain ; 
	tensor2<float,3,3> _commit_strain ; 

	// tensor2<float,3,3>         iterate_stress;                    // Iterative Stress State
	// tensor2<float,3,3>         iterate_strain;                    // Iterative Strain State
	// tensor2<float,3,3>         iterate_plastic_strain;            // Iterative Plastic Strain State

	// tensor2<float,3,3>         converge_commit_stress;            // Commit/Stored Stress State
	// tensor2<float,3,3>         converge_commit_strain;            // Commit/Stored Strain State
	// tensor2<float,3,3>         converge_commit_plastic_strain;    // Commit/Stored Plastic Strain State

	// tensor2<float,3,3>         save_iter_stress;                  // Commit/Stored Stress State
	// tensor2<float,3,3>         save_iter_strain;                  // Commit/Stored Strain State
	// tensor2<float,3,3>         save_iter_plastic_strain;          // Commit/Stored Plastic Strain State




	static tensor4<float,3,3,3,3> _Ee;                    // elastic constant: 3*3*3*3 tensor
	static tensor4<float,3,3,3,3> _D_inv;                 // elastic constant: 3*3*3*3 tensor
	static tensor4<float,3,3,3,3> _Eep;                    // elastic constant: 3*3*3*3 tensor

	static const tensor2<float,3,3> kronecker_delta ;// Delta 

	eindex < 'i' > I;                       // Dummy or Free eindex for LTensor
	eindex < 'j' > J;                       // Dummy or Free eindex for LTensor
	eindex < 'k' > K;                       // Dummy or Free eindex for LTensor
	eindex < 'l' > L;                       // Dummy or Free eindex for LTensor
	eindex < 'o' > O;                       // Dummy or Free eindex for LTensor
	eindex < 't' > T;                       // Dummy or Free eindex for LTensor
	eindex < 'x' > X;                       // Dummy or Free eindex for LTensor
	eindex < 'y' > Y;                       // Dummy or Free eindex for LTensor

	static const  tensor1<float,3> _ZeroVec3;
	static const  tensor2<float,3,3> _ZeroMat3;
	static const  tensor4<float,3,3,3,3> _ZeroStiff;


};