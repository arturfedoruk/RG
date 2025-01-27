#include "smdr.h"
//#include "/home/gebrem/Downloads/SMDR-main/src/smdr_pdg_2021.h"
#include "iostream"
#include "fstream"
#include "string"

#define ZEROSAFE(a) (((a) > (SMDR_TOL)) ? (a) : (SMDR_TOL)) //idk wht's that
using namespace std; //based

//-----------------------------------------------------------------------------------------

//my version of SMDR_Fit_Inputs, where the loops for {Mt, Mt, Mh, MZ, MW, GFermi, QCDQED_at_MZ, mbmb} stored in loop_config
int  my_Fit_Inputs (SMDR_REAL Q_target,
                      SMDR_REAL alphaS_5_MZ_target,
                      SMDR_REAL alpha_target,
                      SMDR_REAL GFermi_target,
                      SMDR_REAL MZ_target,
                      SMDR_REAL Mh_target,
                      SMDR_REAL Mt_target,
                      SMDR_REAL mbmb_target,
                       // i totally eliminated all leptons and quarks lighter than b from the whole code (they're too smol to account for)
                      SMDR_REAL Delta_alpha_had_target,
                      SMDR_REAL error_target,
                      float* loop_config)
{
  char funcname[] = "my_Fit_Inputs";
  int j;
  SMDR_REAL vratio, v2ratio;
  SMDR_REAL KZ, Ka;
  SMDR_REAL err_total, err_alphaS, err_alpha, err_MZ, err_GFermi, err_Mh; 
  SMDR_REAL err_Mt, err_mbmb;  
  SMDR_REAL g3rat, lambdarat, ytrat, ybrat;
  SMDR_REAL del_v, del_g3, del_g, del_gp, del_yt, del_yb, del_lambda;
  int result_niters;
  
  // /home/gebrem/Downloads/ is the directory where your SMDR_main is located. I know that that's not very demure and is kinda crutchy
  #include "/home/gebrem/Downloads/SMDR-main/src/accelcoeffs.h"

  if (mbmb_target < 3.0) {
    SMDR_Warn (funcname, "Target mb(mb) is too small.");
    SMDR_Warn (funcname, "Setting it and all lighter fermion masses to 0.\n");
    mbmb_target = 0;
  } 

  /* Don't compute Mt, Mh with needlessly high tolerances... */ 
  SMDR_MTPOLE_TOLERANCE = 1.0e-8;
  SMDR_MHPOLE_TOLERANCE = 1.0e-8;
  if (error_target < 1.0e-7) {
    SMDR_MTPOLE_TOLERANCE = error_target/10.;
    SMDR_MHPOLE_TOLERANCE = error_target/10.;
  }
  if (SMDR_MTPOLE_TOLERANCE < 1.01e-15) SMDR_MTPOLE_TOLERANCE = 1.01e-15;
  if (SMDR_MHPOLE_TOLERANCE < 1.01e-15) SMDR_MHPOLE_TOLERANCE = 1.01e-15;

  /* Use ReferenceModel.dat as our initial guess for the MSbar
     parameters, and their corresponding OS predictions.
  */
  /* SMDR_Read_MSbar_Inputs ("ReferenceModel.dat"); */
  /* SMDR_Read_Value ("ReferenceModel.dat", "SMDR_v_in"); */
  /* SMDR_Read_OS_Inputs ("ReferenceModel.dat"); */
  #include "/home/gebrem/Downloads/SMDR-main/src/AUTOrefmodel.h"

  SMDR_RGeval_SM (Mt_target, loop_config[6]);

  SMDR_Save_Inputs ();
  SMDR_Delta_alpha_had_5_MZ_in = Delta_alpha_had_target;

  for (j=0; j<18; j++) {

    /* Find the new input parameter estimates by comparing OS observables to the 
       targets. On the 0th iteration, the OS observables and MSbar parameters
       are both as specified in "ReferenceModel.dat".
    */
    v2ratio = SMDR_GFermi/GFermi_target;
    vratio = SMDR_SQRT (v2ratio);

    KZ = (SMDR_g_in * SMDR_g_in + SMDR_gp_in * SMDR_gp_in) * 
         (MZ_target * MZ_target)/
         (SMDR_MZ_PDG * SMDR_MZ_PDG * v2ratio);

    Ka = ((SMDR_g_in * SMDR_g_in * SMDR_gp_in * SMDR_gp_in)/
          (SMDR_g_in * SMDR_g_in + SMDR_gp_in * SMDR_gp_in)) * 
          (alpha_target/SMDR_alpha);
	
	//here I made g_in -> SMDR_g_in and gp_in -> SMDR_gp_in
    del_g  = (SMDR_SQRT(KZ * (1. + SMDR_SQRT(1. - 4. * Ka/KZ))/2.)/SMDR_g_in) - 1.;
    del_gp = (SMDR_SQRT(KZ * (1. - SMDR_SQRT(1. - 4. * Ka/KZ))/2.)/SMDR_gp_in) - 1.;

    del_g3 = SMDR_SQRT (alphaS_5_MZ_target/SMDR_alphaS_5_MZ) - 1.;
    del_v = vratio - 1.;
    del_lambda = ((Mh_target * Mh_target)/
                 (SMDR_Mh_pole * SMDR_Mh_pole * v2ratio)) - 1.;
    del_yt = (Mt_target/(SMDR_Mt_pole * vratio)) - 1.;
    del_yb = (mbmb_target/(ZEROSAFE(SMDR_mbmb) * vratio)) - 1.;

    g3rat = 1. + cog3_g3 * del_g3 
               + cog3_g32 * del_g3 * del_g3 
               + cog3_yt * del_yt;

    lambdarat = 1. + del_lambda 
                   + colambda_yt * del_yt 
                   + colambda_g3 * del_g3;

    ytrat = 1. + coyt_yt * del_yt 
               + coyt_yt2 * del_yt * del_yt
               + coyt_g3 * del_g3;

    ybrat = 1. + coyb_yb * del_yb 
               + coyb_yb2 * del_yb * del_yb +
               + coyb_g3 * del_g3
               + coyb_yt * del_yt;

    if ((ybrat < 0) || (mbmb_target < SMDR_TOL)) ybrat = 0;

    SMDR_v_in *= vratio;
    SMDR_g3_in *= g3rat;
    SMDR_g_in  *= 1. + del_g;
    SMDR_gp_in *= 1. + del_gp;
    SMDR_lambda_in *= lambdarat;
    SMDR_yt_in *= ytrat;
    SMDR_yb_in *= ybrat;

    err_alphaS = SMDR_alphaS_5_MZ/alphaS_5_MZ_target - 1.;
    err_alpha = SMDR_alpha/alpha_target - 1.;
    err_MZ = SMDR_MZ_PDG/MZ_target - 1.;
    err_GFermi = SMDR_GFermi/GFermi_target - 1.;
    err_Mh = SMDR_Mh_pole/Mh_target - 1.;
    err_Mt = SMDR_Mt_pole/Mt_target - 1.;
    err_mbmb = (SMDR_mbmb - mbmb_target)/SMDR_mbmb_EXPT;

    err_total = SMDR_FABS(err_alphaS) +
                SMDR_FABS(err_alpha) + 
                SMDR_FABS(err_MZ) + 
                SMDR_FABS(err_GFermi) + 
                SMDR_FABS(err_Mh) + 
                SMDR_FABS(err_Mt) + 
                SMDR_FABS(err_mbmb);
                
      // i deleted some commented sectors for debugging, for full version see fit_obc.c

    if (j>0) {
      printf("Iteration %d: ",j);
      printf("Total fractional error = %.15Lf\n", err_total); 
      if (err_total < error_target) break;
    }
  
    /* Now evaluate all of the on-shell observables. */
    SMDR_Eval_Mt_pole (173.22, 1, loop_config[0], loop_config[1], &SMDR_Mt_pole, &SMDR_Gammat_pole);
    SMDR_Eval_Mh_pole (173.22, loop_config[2], &SMDR_Mh_pole, &SMDR_Gammah_pole);
    SMDR_Eval_MZ_pole (173.22, loop_config[3], &SMDR_MZ_pole, &SMDR_GammaZ_pole,
                       &SMDR_MZ_PDG, &SMDR_GammaZ_PDG);

    /* This one is computed only because it is needed by SMDR_Eval_Gauge() */
    SMDR_Eval_MW_pole (173.22, loop_config[4], &SMDR_MW_pole, &SMDR_GammaW_pole,
                       &SMDR_MW_PDG, &SMDR_GammaW_PDG);

    SMDR_GFermi = SMDR_Eval_GFermi (SMDR_Mt_EXPT, loop_config[5]);
    SMDR_Eval_Gauge (SMDR_Mt_pole, SMDR_Mh_pole, SMDR_MW_PDG);
    SMDR_Eval_QCDQED_at_MZ (SMDR_MZ_PDG, 173.22, loop_config[6]);

    if (mbmb_target > SMDR_TOL) 
      SMDR_mbmb = SMDR_Eval_mbmb(173.22, loop_config[7]);
    else SMDR_mbmb = 0;

    result_niters = j;
  }

  SMDR_Load_Inputs();

  SMDR_Lambda = 0;
  SMDR_Lambda_in = SMDR_Lambda = -SMDR_Eval_Veffmin (-1, loop_config[2]);
  SMDR_m2_in = SMDR_m2;

  SMDR_RGeval_SM (Q_target, loop_config[6]);
  SMDR_Save_Inputs();

  return(result_niters);
};


//----------------------------------------------------------------------------------------


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
