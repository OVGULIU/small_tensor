#include "vonMisesLinearHardening.h"


namespace openus{

vonMisesLinearHardening::vonMisesLinearHardening(){

}


vonMises::Mat33 
vonMisesLinearHardening::kinematic_derivative(Mat33 const& m, Mat33 const& back_stress) const {
	Mat33 sol; 
	float mm = - 1./3. * ( m(0,0) + m(1,1) + m(2,2) ) ;
	sol(I,J) = _kinematic_harden_rate * ( m(I,J) + mm * kronecker_delta(I,J) ) ;
	return sol ;
}


} // namespace openus

