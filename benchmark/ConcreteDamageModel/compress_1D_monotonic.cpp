#include "ConcreteDamageModel.h"
#include <fstream>

using namespace smalltensor; 
int main(int argc, char const *argv[])
{
	// double bulk_modulus = (double) atof(argv[1]);
	// double scale_hardening = (double) atof(argv[2]);
	// double max_strain_in = (double) atof(argv[3]);
	// double strain_incr = (double) atof(argv[4]);
	// int Nloop   = atoi(argv[5]);

	ofstream outfile;
	outfile.open("strain_stress.txt");

	int material_tag{1}  ;

	float E_in = 26e9 ; 
	float v_in = 0.16f ; 
	float beta_in = 0.59f ; // plastic deformation rate
	float K_in = 0.49f ; 
	float r_0_tensile_in = 300 ; 
	float r_0_compress_in = 4000 ; 
	float mesh_A_tensile_in = 0.1f ; 
	float mesh_A_compress_in = 1.5f ; 
	float mesh_B_compress_in = 0.75f ; 

	auto theMaterial= new ConcreteDamageModel(
		E_in,
		v_in,
		beta_in,
		K_in,
		r_0_tensile_in,
		r_0_compress_in,
		mesh_A_tensile_in,
		mesh_A_compress_in,
		mesh_B_compress_in
	);

	float max_strain_in = 0.01   ;
	float strain_incr = 1E-5  ;
	auto stress_ret = theMaterial->getStressTensor();
	auto strain_ret = theMaterial->getStrainTensor();

	outfile << strain_ret(0,1) <<"\t" 
			<< stress_ret(0,1) <<"\t"
			<< endl ;

	int Nsteps = max_strain_in/strain_incr ;

	// Loading 
	tensor2<float,3,3> input_strain ;
	for (int i = 0; i < Nsteps; ++i){
		cout << " -------------------------------------------------- " << endl;
		cout << " - step " << i << endl;
		
		// input_strain(0,0) =  incr_size;
		input_strain(0,0) = - strain_incr ;
		// input_strain(1,1) = - strain_incr ;
		// input_strain(2,2) = - strain_incr ;

		theMaterial->setTrialStrainIncr(input_strain);

		theMaterial->commitState();
		stress_ret = theMaterial->getStressTensor();
		strain_ret = theMaterial->getStrainTensor();

		outfile << - strain_ret(0,0)  <<"\t" 
				<< - stress_ret(0,0) <<"\t"
				<< endl ;
	}

	// Done the experiments and Clean
	delete theMaterial ;
	outfile.close();

	return 0;

}











// for (int loop = 0; loop < Nloop; ++loop)
// {
// 	// unloading 
// 	for (int i = 0; i < 2 * Nsteps; ++i)
// 	{
// 		cout<< "--------------------------------------------------" <<endl;
// 		cout<< "unloading step " << i <<endl;
// 		input_strain *= 0. ;
// 		// input_strain(0,0) = - strain_incr;
// 		input_strain(0,1) = - strain_incr / 2.;
// 		input_strain(1,0) = - strain_incr / 2.;

// 		theMaterial->setTrialStrainIncr(input_strain);

// 		theMaterial->commitState();
// 		stress_ret = theMaterial->getStressTensor();
// 		strain_ret = theMaterial->getStrainTensor();
// 		Nactive = theMaterial->getNumActiveYS();

// 		outfile << strain_ret(0,1)  <<"\t" 
// 				<< stress_ret(0,1) <<"\t"
// 				<< Nactive << endl ;
// 	}

// 	// reloading 
// 	for (int i = 0; i < 2 * Nsteps; ++i)
// 	{
// 		cout<< "--------------------------------------------------" <<endl;
// 		cout<< "reloading step " << i <<endl;
// 		input_strain *= 0. ;
// 		// input_strain(0,0) =  strain_incr;
// 		input_strain(0,1) = strain_incr / 2.;
// 		input_strain(1,0) = strain_incr / 2.;

// 		theMaterial->setTrialStrainIncr(input_strain);

// 		theMaterial->commitState();
// 		stress_ret = theMaterial->getStressTensor();
// 		strain_ret = theMaterial->getStrainTensor();

// 		outfile << strain_ret(0,1)  <<"\t" 
// 				<< stress_ret(0,1) <<"\t"
// 				<< endl ;
// 	}
// }