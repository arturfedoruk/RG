#include "smdr.h"
#include "iostream"
#include "fstream"
#include "string"

using namespace std;

#include "my_Fit_Inputs.cpp"
#include "loop_configs.cpp"

int main(){

	#include "smdr_pdg_2025.h"	
	
	ofstream output_file("for_theorerr.txt");
	
	double Q0 = 173.22;
	
	double ERROR_TOLERANCE = 1.e-6; // idk
	double pi = 3.14159265359;
	
	float factors[5] = {1/2., 1/1.5, 1., 1.5, 2.};
	
	for(int iconfig = 0; iconfig < nconfigs; iconfig++){
	
		output_file << "\n++++++++++++++++++++++++++++++++++++++++\n" << loop_names[iconfig] << endl; 
		cout << loop_names[iconfig] << endl;
	
		for(int ifactor = 0; ifactor < 5; ifactor++){
		
			SMDR_Q_in = Q0 * factors[ifactor];
			
			#include "smdr_pdg_2025.h"
			
			my_Fit_Inputs (SMDR_Q_in,
				   SMDR_alphaS_MZ_EXPT,
				   SMDR_alpha_EXPT,
				   SMDR_GFermi_EXPT,
				   SMDR_MZ_EXPT,
				   SMDR_Mh_EXPT,
				   SMDR_Mt_EXPT,
				   SMDR_mbmb_EXPT,
				   SMDR_Delta_alpha_had_5_MZ_EXPT,
				   ERROR_TOLERANCE,
				   loop_configs[iconfig]);
				   
			 output_file << "Q = " << Q0 * factors[ifactor] << "\n x1 = " << SMDR_gp_in*SMDR_gp_in*5 / (48*pi*pi) << ", x2 = " << SMDR_g_in*SMDR_g_in / (16*pi*pi) << ", x3 = " << SMDR_g3_in*SMDR_g3_in / (16*pi*pi) << ", \n y1 = " << SMDR_yt_in*SMDR_yt_in / (16*pi*pi) << ", y2 = " << SMDR_yb_in*SMDR_yb_in / (16*pi*pi) << ", z = " << SMDR_lambda_in / (16*pi*pi) << endl;
			                    
                   }
                   
       }
	
	return 0;
	
}
