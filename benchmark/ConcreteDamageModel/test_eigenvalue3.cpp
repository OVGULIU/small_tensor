#include "ConcreteDamageModel_utils.hpp"

#include <iostream>

using namespace smalltensor ; 

int main(int argc, char const *argv[])
{
	tensor2<float,3,3> matrix ; 
	matrix(0,0) = 1 ; matrix(0,1) = 2 ; matrix(0,2) = 3 ; 
	matrix(1,0) = 2 ; matrix(1,1) = 4 ; matrix(1,2) = 5 ; 
	matrix(2,0) = 3 ; matrix(2,1) = 5 ; matrix(2,2) = 6 ; 

	tensor1<float,3> eigenvalues ;

	eigenvalues = getEigenValues3(matrix) ; 

	std::cerr << " eigenvalues = "<< eigenvalues << std::endl ;

	return 0;
}