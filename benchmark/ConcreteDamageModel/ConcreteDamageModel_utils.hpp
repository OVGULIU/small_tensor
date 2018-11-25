
#include "../../smalltensor/smalltensor.h"
using namespace smalltensor;

#include <cmath>
#include <iostream>

// Input: 3x3 real-symmetric matrix
// Output: 3x1 eigenvalues eig3 <= eig2 << eig1

// Reference: https://en.wikipedia.org/wiki/Eigenvalue_algorithm

tensor1<float,3> GetEigenValues3(tensor2<float, 3,3> const& matrix){
	eindex < 'i' > I;                       // Dummy or Free eindex for LTensor
	eindex < 'j' > J;                       // Dummy or Free eindex for LTensor
	const tensor2<float,3,3> kronecker_delta("identity");

	tensor1<float,3> result ; 

	const float trace = matrix(0,0) + matrix(1,1) + matrix(2,2) ;
	const float p1 = std::pow(matrix(0,1),2) + std::pow(matrix(0,2),2) + std::pow(matrix(1,2),2) ; 

	if ( std::fabs(p1) < 1e-12 ){
		// diagonal matrix 
		result(0) = std::max( std::max(matrix(0,0), matrix(1,1)), matrix(2,2) ) ; // max-eigenvalue
		result(2) = std::min( std::min(matrix(0,0), matrix(1,1)), matrix(2,2) ) ; // min-eigenvalue
		result(1) = trace - result(0) - result(2) ; // mid-eigenvalue
	}else{
		const float q = trace /3.f ; 
		const float p2 = std::pow( matrix(0,0) - q, 2) + std::pow( matrix(1,1) - q, 2) 
					   + std::pow( matrix(2,2) - q, 2) + 2 * p1 ; 
		const float p = std::sqrt( p2 / 6.f ) ;
		tensor2<float,3,3> B ; 
		B(I,J) = (1.f / p ) * ( matrix(I,J) - q * kronecker_delta(I,J) ) ; 
		const float r = B.compute_Determinant() / 2.f ; 

		// In exact arithmetic for a symmetric matrix -1 <= r <= 1
		// but computation error can leave it slightly outside this range. 
		float phi = 0.f ; 
		if ( r <= -1){
			phi = M_PI / 3.f ; 
		}else if ( r >= 1 ){
			phi = 0 ; 
		}else{
			phi = std::acos(r) / 3.f ;
		}

		// the eigenvalues satisfy eig3 < eig2 < eig1 
		result(0) = q + 2 * p * std::cos(phi) ; 
		result(2) = q + 2 * p * cos( phi + (2 * M_PI/ 3) ) ; 
		result(1) = 3 * q - result(0) - result(2) ; 

	}
	return result ;
}




// 
float Macaulay(float value){
	return value > 0 ? value : 0 ; 
}

float Heaviside(float value){
	return value >= 0 ? 1 : 0 ; 
}



void SplitStress(tensor2<float,3,3> const& input, 
	tensor2<float,3,3>& tensile_part, 
	tensor2<float,3,3>& compress_part, 
	)
{

	tensor1<float,3> eigenvalues = GetEigenValues3(input) ; 

	tensile_part.zero() ;
	tensile_part(0,0) = Macaulay(eigenvalues(0)) ; 
	tensile_part(1,1) = Macaulay(eigenvalues(1)) ; 
	tensile_part(2,2) = Macaulay(eigenvalues(2)) ; 

	compress_part(I,J) = input(I,J) - tensile_part(I,J) ; 

}


float GetEquivStressCompress(tensor2<float,3,3> const& input, float K_compress){
	tensor1<float,3> eigenvalues = GetEigenValues3(input) ; 

	const float octahedral_nomral_stress = 1.f/3.f * ( eigenvalues(0) + eigenvalues(1) + eigenvalues(2) ) ; 
	const float octahedral_shear_stress = 1.f/3.f * std::sqrt(
					std::pow( eigenvalues(0) - eigenvalues(1), 2)
				  + std::pow( eigenvalues(1) - eigenvalues(2), 2)
				  + std::pow( eigenvalues(0) - eigenvalues(2), 2)
				);

	const float compress_equiv_stress = std::sqrt( std::sqrt(3) * 
			(K_compress * octahedral_nomral_stress + octahedral_shear_stress)
		) ; 

	return compress_equiv_stress ;
}