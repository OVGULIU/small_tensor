#include "vonMisesArmstrongFrederick.h"


namespace openus{

vonMisesArmstrongFrederick::vonMisesArmstrongFrederick(){

}


vonMisesArmstrongFrederick::Mat33 
vonMisesArmstrongFrederick::kinematic_derivative(Mat33 const& m, Mat33 const& back_stress) const {
	Mat33 mdev; 
	float mm = - 1./3. * ( m(0,0) + m(1,1) + m(2,2) ) ;
	mdev(I,J) = m(I,J) - mm * kronecker_delta(I,J) ;

	Mat33 sol; 
	sol(I,J) = 2./3 * _AF_ha * mdev(I,J) - _AF_cr * sqrt( 2./3. * mdev(K,L) * mdev(K,L) ) * back_stress(I,J) ;

	return sol ;
}


} // namespace openus

