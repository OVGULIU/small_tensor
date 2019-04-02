#include "vonMises.h"
#include <fstream>
#include <memory>

int main(int argc, char const *argv[])
{
	using namespace openus;
	using namespace smalltensor ;
	using namespace std;
	ofstream outfile;
	outfile.open("strain_stress.txt");

	auto material = std::make_shared<vonMises>() ;
	material->Initialize();

	tensor2<float,3,3> input_strain ;

	input_strain *= 0.;
	float strain_incr = 0.0005 ; 
	input_strain(0,1) = strain_incr;
	input_strain(1,0) = strain_incr;

	constexpr int Nsteps = 20 ;

	std::cout << material->_iter_strain(0,1) << " " << material->_iter_stress(0,1) << std::endl; 
	outfile << material->_iter_strain(0,1) << " " << material->_iter_stress(0,1) << std::endl; 
	for (int i = 0; i < Nsteps; ++i)	{
		material->SetTrialStrainIncr(input_strain) ;
		material->CommitState();
		std::cout << material->_iter_strain(0,1) << " " << material->_iter_stress(0,1) << std::endl; 
		outfile << material->_iter_strain(0,1) << " " << material->_iter_stress(0,1) << std::endl; 
	}

	outfile.close();
	return 0;
}


