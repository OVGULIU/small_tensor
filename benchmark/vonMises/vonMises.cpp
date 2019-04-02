#include "vonMises.h"


namespace openus{

const vonMises::Mat33 vonMises::_zero_strain ;
const vonMises::Mat33 vonMises::_zero_stress ;
const vonMises::Mat33 vonMises::kronecker_delta("identity") ;

vonMises::vonMises()
{

}

void 
vonMises::Initialize() {
	_stiff_tensor *= 0. ;
	auto mu = _E / ( 2 * ( 1 + _nu ) );
	auto lambda = ( _nu * _E ) / ( ( 1 + _nu ) * ( 1 - 2 * _nu ) );
	_stiff_tensor( 0, 0, 0, 0 ) = lambda + 2 * mu;
	_stiff_tensor( 0, 0, 1, 1 ) = lambda;
	_stiff_tensor( 0, 0, 2, 2 ) = lambda;
	_stiff_tensor( 0, 1, 0, 1 ) = mu;
	_stiff_tensor( 0, 1, 1, 0 ) = mu;
	_stiff_tensor( 0, 2, 0, 2 ) = mu;
	_stiff_tensor( 0, 2, 2, 0 ) = mu;
	_stiff_tensor( 1, 0, 0, 1 ) = mu;
	_stiff_tensor( 1, 0, 1, 0 ) = mu;
	_stiff_tensor( 1, 1, 0, 0 ) = lambda;
	_stiff_tensor( 1, 1, 1, 1 ) = lambda + 2 * mu;
	_stiff_tensor( 1, 1, 2, 2 ) = lambda;
	_stiff_tensor( 1, 2, 1, 2 ) = mu;
	_stiff_tensor( 1, 2, 2, 1 ) = mu;
	_stiff_tensor( 2, 0, 0, 2 ) = mu;
	_stiff_tensor( 2, 0, 2, 0 ) = mu;
	_stiff_tensor( 2, 1, 1, 2 ) = mu;
	_stiff_tensor( 2, 1, 2, 1 ) = mu;
	_stiff_tensor( 2, 2, 0, 0 ) = lambda;
	_stiff_tensor( 2, 2, 1, 1 ) = lambda;
	_stiff_tensor( 2, 2, 2, 2 ) = lambda + 2 * mu;
}

void 
vonMises::Update() {

}

vonMises::Tensor4 const&
vonMises::GetStiffnessTensor() const {
	return _stiff_tensor ;
}

int 
vonMises::SetTrialStrainIncr(Mat33 const& strain_incr){
	_iter_strain(I,J) = _commit_strain(I,J) + strain_incr(I,J) ;
	return compute_stress(strain_incr) ;
}

int
vonMises::compute_stress(Mat33 const& strain_incr){
	Mat33 stress_incr;
	stress_incr(I,J) = _stiff_tensor(I,J,K,L) * strain_incr(K,L) ;

	Mat33 predict_stress ;
	predict_stress(I,J) = _commit_stress(I,J) + stress_incr(I,J) ;

	float yf_val_start = yield_surface_val(_commit_stress, _iter_back_stress, _iter_yf_radius ) ; 
	float yf_val_end = yield_surface_val(predict_stress, _iter_back_stress, _iter_yf_radius ) ; 
	_iter_stress(I,J) = predict_stress(I,J) ;

	if ( (yf_val_start <= 0.0 && yf_val_end <= 0.0) || yf_val_start > yf_val_end ){
		// elastic
		return 0 ;
	}else{
		// plastic
		auto n = df_dsigma(_iter_stress, _iter_back_stress) ;
		auto m = n ;
		float hardening = hardening_ksi_h( _iter_stress , _iter_back_stress, m ) ;
		float denominator = n(I,J) * _stiff_tensor(I,J,K,L) * m(K,L) - hardening ;
		float dLambda = yf_val_end /denominator ;
		_iter_stress(I,J) = predict_stress(I,J) - dLambda * _stiff_tensor(I,J,K,L) * m(K,L) ;
		// internal variables evolve.
		evolve_internal_variables(dLambda, _iter_back_stress, _iter_yf_radius, m) ;
	}
	return 0;
}

float 
vonMises::yield_surface_val(Mat33 const& stress, Mat33 const& back_stress, float radius) const {
	float sol = 0 ;
	float pp = -1./3 * (stress(0,0) + stress(1,1) + stress(2,2) ) ;
	Mat33 s ;
	s(I,J) = stress(I,J) + pp * kronecker_delta(I,J) ;
	Mat33 s_minus_alpha ;
	s_minus_alpha(I,J) = s(I,J) - back_stress(I,J) ;
	float tmp = sqrt( s_minus_alpha(I,J) * s_minus_alpha(I,J) ) ; 
	sol = tmp - sqrt(2./3.) * radius ;
	return sol ;
}

vonMises::Mat33 
vonMises::df_dsigma(Mat33 const& stress, Mat33 const& back_stress) const {
	float pp = -1./3 * (stress(0,0) + stress(1,1) + stress(2,2) ) ;
	Mat33 s ;
	s(I,J) = stress(I,J) + pp * kronecker_delta(I,J) ;
	Mat33 s_minus_alpha ;
	s_minus_alpha(I,J) = s(I,J) - back_stress(I,J) ;
	float denominator = sqrt( s_minus_alpha(I,J) * s_minus_alpha(I,J) ) ; 
	s_minus_alpha(I,J) = s_minus_alpha(I,J) / denominator ;
	return s_minus_alpha ;
}

float 
vonMises::hardening_ksi_h(Mat33 const& stress, Mat33 const& back_stress, Mat33 const& m) const {
	float pp = -1./3 * (stress(0,0) + stress(1,1) + stress(2,2) ) ;
	Mat33 s ;
	s(I,J) = stress(I,J) + pp * kronecker_delta(I,J) ;
	Mat33 s_minus_alpha ;
	s_minus_alpha(I,J) = s(I,J) - back_stress(I,J) ;
	float denominator = sqrt( s_minus_alpha(I,J) * s_minus_alpha(I,J) ) ; 
	if ( std::fabs(denominator) < 1e-5 ){
		return 0 ;
	}
	// isotropic
	float sol = 0 ;
	sol = - sqrt(2./3.) * isotropic_derivative(m) ;

	// kinematic
	Mat33 alpha_deriv = kinematic_derivative(m) ;
	sol = - s_minus_alpha(I,J) / denominator * kinematic_derivative(m)(I,J) ; 
	return sol ;
}

float 
vonMises::isotropic_derivative(Mat33 const& m) const {
	float m_eq = sqrt( 2./3. * m(I,J) * m(I,J) ) ;
	return _isotropic_harden_rate * m_eq ;
}

vonMises::Mat33 
vonMises::kinematic_derivative(Mat33 const& m) const {
	Mat33 sol; 
	float mm = - 1./3. * ( m(0,0) + m(1,1) + m(2,2) ) ;
	sol(I,J) = _kinematic_harden_rate * ( m(I,J) + mm * kronecker_delta(I,J) ) ;
	return sol ;
}

void
vonMises::evolve_internal_variables(float dLambda, 
	Mat33& _iter_back_stress, float & _iter_yf_radius, Mat33 const& m) const {
	_iter_back_stress(I,J) = _iter_back_stress(I,J) + dLambda * kinematic_derivative(m)(I,J) ;
	_iter_yf_radius += dLambda * isotropic_derivative(m) ;
}

void 
vonMises::CommitState(){
	_commit_stress         = _iter_stress ;
	_commit_strain         = _iter_strain ;
	_commit_plastic_strain = _iter_plastic_strain ;
	_commit_back_stress    = _iter_back_stress ;
}


} // namespace openus

