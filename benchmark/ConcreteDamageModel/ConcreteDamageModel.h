
#include "../../smalltensor/smalltensor.h"
#include <limits>
#include <vector>
#include <cmath>

using namespace smalltensor ; 
class ConcreteDamageModel
{
public:
	ConcreteDamageModel(

		);
	// ~ConcreteDamageModel();
	

	int setTrialStrainIncr(tensor2<float,3,3> const& strain_incr);
	void commitState();
	tensor2<float,3,3> getStressTensor() const;
	tensor2<float,3,3> getStrainTensor() const;


private:
	int compute_effective_stress() ; 
	int compute_cauchy_stress() ; 
private:
	float _d_tensile ; 
	float _d_compress ; 

	float _r_0_tensile ; 
	float _r_0_compress ; 
	
	float _r_tensile ; 
	float _r_compress ; 

	float _E;                              // Elastic Modulus
	float _v;                              // Poisson Ratio
	float _rho;                            // Density

	tensor2<float,3,3> _effective_stress ;                    // 

	tensor2<float,3,3> _effective_stress_trial ;
	tensor2<float,3,3> _effective_stress_commit ;


	tensor2<float,3,3> _effective_stress_compress;            // -
	tensor2<float,3,3> _effective_stress_tensile;             // +



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