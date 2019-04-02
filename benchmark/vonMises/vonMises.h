#ifndef OPENUS_FE_MATERIAL_vonMises_H_
#define OPENUS_FE_MATERIAL_vonMises_H_

#include <smalltensor/smalltensor.h>
#include <cmath>

namespace openus{

class vonMises
{
typedef smalltensor::tensor2<float,3,3> Mat33 ; 
typedef smalltensor::tensor4<float,3,3,3,3> Tensor4 ;

public:
	vonMises();
	void Initialize() ; 
	void Update() ; 

	int SetTrialStrainIncr(Mat33 const& strain_incr) ;
	Tensor4 const& GetStiffnessTensor() const ;
	void CommitState();

	int compute_stress(Mat33 const& strain_incr) ;
	float yield_surface_val(Mat33 const& stress, Mat33 const& back_stress, float radius) const ;
	Mat33 df_dsigma(Mat33 const& stress, Mat33 const& back_stress) const ;
	float hardening_ksi_h(Mat33 const& stress, Mat33 const& back_stress, Mat33 const& m) const ;
	float isotropic_derivative(Mat33 const& m) const ;
	Mat33 kinematic_derivative(Mat33 const& m) const ;

	float _E = 1e3 ;
	float _nu = 0.0 ;
	float _isotropic_harden_rate = 3 ;
	float _kinematic_harden_rate = 4 ;
	float _iter_yf_radius = 2 ; 
	float _commit_yf_radius = 2 ; 

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

