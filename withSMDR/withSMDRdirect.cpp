#include "smdr.h"
#include "iostream"
#include "fstream"
#include "string"

#define ZEROSAFE(a) (((a) > (SMDR_TOL)) ? (a) : (SMDR_TOL)) //idk wht's that
using namespace std; //based

#include "my_Fit_Inputs.cpp"


int main(){

	ofstream output_file("withSMDRdirectOutput.txt");
	output_file << "\n++++++++++++++++++++++++++++++++++++++++\n DIRECT ALGORITHM \n\n";
	
	SMDR_REAL temp1 = 0;
	SMDR_REAL temp2 = 0; // idc
	double ERROR_TOLERANCE = 1.e-12; // idk
	double pi = 3.14159265359;
	
	// loop configurations for the functions below: 
	// {Mt, Mt, Mh, MZ, MW, GFermi, QCDQED_at_MZ, mbmb} 
	// {y, ?, z, x, x, z, x, y} -1 loop
	// (Mt has two loop configurations)
	float config_111111[9] = {0, 0, 0, 0, 0, 0, 1, 1} ; // for QCDQED_at_MZ & mbmb loop 0 doesn't exist
	float config_222222[9] = {1, 1, 1, 1, 1, 1, 1, 1} ;
	float config_333333[9] = {2, 2, 2, 2, 2, 2, 2, 2} ;
	float config_333221[9] = {2, 2, 0, 2, 2, 0, 2, 1} ;
	float config_444332[9] = {3, 2, 1, 2.5, 2.5, 1, 3, 2} ; 
	float config_444333[9] = {3, 2, 2, 2.5, 2.5, 2, 3, 2} ; // hfor MZ & MW loop 3 doesn't exist
	
	const int nconfigs = 6;
	float* loop_configs[nconfigs] = {config_111111, config_222222, config_333333, config_333221, config_444332, config_444333};
	string loop_names[nconfigs] = {"(1,1,1;1,1;1)","(2,2,2;2,2;2)","(3,3,3;3,3;3)","(3,3,3;2,2;1)","(4,4,4;3,3;2)","(4,4,4;3,3;3)"};
	
  	SMDR_Start_Timer();
	
	for(int iconfig = 0; iconfig < nconfigs; iconfig++){
	
		SMDR_Q_in = 173.22; // mass of top-quark, we work here
		
		SMDR_REAL SMDR_alphaS_5_MZ_EXPT = 0.1179; // SMDR suggests this.... but what the incertainty then??
		SMDR_alpha_EXPT = 0.0072973525693;
		SMDR_GFermi_EXPT = 0.000011663787;
		SMDR_MZ_EXPT = 91.1876; 
		SMDR_Mh_EXPT = 125.25;
		SMDR_Mt_EXPT = 172.5;
		SMDR_mbmb_EXPT = 4.18; // from smdr_pdg_2021.h, there from Review of Particle Properties 2020 
		
		// my_Fit_Inputs allows loop customization and is quicker due zero mass of unnecessary particles
		my_Fit_Inputs (SMDR_Q_in,
		           SMDR_alphaS_5_MZ_EXPT,
		           SMDR_alpha_EXPT,
		           SMDR_GFermi_EXPT,
		           SMDR_MZ_EXPT,
		           SMDR_Mh_EXPT,
		           SMDR_Mt_EXPT,
		           SMDR_mbmb_EXPT,
		           SMDR_Delta_alpha_had_5_MZ_EXPT,
		           ERROR_TOLERANCE,
		           loop_configs[iconfig]);
		           
	         cout << "\n++++++++++++++++++++++++++++++++++++++++\nComputation of MSbar variables for loop configuration " << loop_names[iconfig] << " at Q=Mt resulted in\n g = " << SMDR_gp_in << ", g' = " << SMDR_g_in << ", g3 = " << SMDR_g3_in << ", \n yt = " << SMDR_yt_in << ", yb = " << SMDR_yb_in << ", lambda = " << SMDR_lambda_in << "\n, which corresponds for values at t=0 of \n x1 = " << SMDR_gp_in*SMDR_gp_in*5 / (48*pi*pi) << ", x2 = " << SMDR_g_in*SMDR_g_in / (16*pi*pi) << ", x3 = " << SMDR_g3_in*SMDR_g3_in / (16*pi*pi) << ", \n y1 = " << SMDR_yt_in*SMDR_yt_in / (16*pi*pi) << ", y2 = " << SMDR_yb_in*SMDR_yb_in / (16*pi*pi) << ", z = " << SMDR_lambda_in / (16*pi*pi) << "\n++++++++++++++++++++++++++++++++++++++++\n" << endl;
	         
        	output_file << loop_names[iconfig] << "\n x1 = " << SMDR_gp_in*SMDR_gp_in*5 / (48*pi*pi) << ", x2 = " << SMDR_g_in*SMDR_g_in / (16*pi*pi) << ", x3 = " << SMDR_g3_in*SMDR_g3_in / (16*pi*pi) << ", \n y1 = " << SMDR_yt_in*SMDR_yt_in / (16*pi*pi) << ", y2 = " << SMDR_yb_in*SMDR_yb_in / (16*pi*pi) << ", z = " << SMDR_lambda_in / (16*pi*pi) << "\n";
                   
       }
       
       output_file << "\n++++++++++++++++++++++++++++++++++++++++\n";
       
       SMDR_Timer();
       cout  << "(what's the matter with g and g'? why no consistency between Robertson and Bednyakov?)\n\n" << "Total calculation time: " << SMDR_Time_Total << " sec" << endl;

	return 0;
}
