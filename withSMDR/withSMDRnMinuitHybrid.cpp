// to launch the program: 
// g++ `root-config --cflags` withSMDRnMinuitHybrid.cpp `root-config --libs` -lm -lsmdr -ltsil -l3vil

#include "smdr.h"
#include "iostream"
#include "fstream"
#include "string"
#include "TROOT.h"
#include "TF1.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TError.h"
using namespace std;
#define ZEROSAFE(a) (((a) > (SMDR_TOL)) ? (a) : (SMDR_TOL)) //idk wht's that

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
      /*printf("Iteration %d: ",j);
      printf("Total fractional error = %.15Lf\n", err_total); */
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


// function, that returns specific MSbar quantity
Double_t compute_MSbar  (const Double_t* pars) // Q_target, alphaS_5_MZ_target, alpha_target, GFermi_target, MZ_target, Mh_target, Mt_target, mbmb_target, Delta_alpha_had_target, error_target, loop_config0, loop_config1, loop_config2, loop_config3, loop_config4, loop_config5, loop_config6, loop_config7, what (what you want to compute)
{
	
	float loop_config[8] = {pars[10],pars[11],pars[12],pars[13],pars[14],pars[15],pars[16],pars[17]};
        my_Fit_Inputs (pars[0], pars[1], pars[2], pars[3], pars[4], pars[5], pars[6], pars[7], pars[8], pars[9], loop_config);        
                      
	switch (int(pars[18])){
	case 1 : return SMDR_gp_in; break; // g
	case 2 : return SMDR_g_in; break; // g'
	case 3 : return SMDR_g3_in; break; // g3
	case 4 : return SMDR_yt_in; break; // yt
	case 5 : return SMDR_yb_in; break; // yb
	case 6 : return SMDR_lambda_in; break; // lambda
	default: cout << "!!!invalid 'what' parameter!!!" << endl; return -999;
	}                   
}

// the same but with minus (needed for maximization)
Double_t compute_MSbar_butminus  (const Double_t* pars)
{
	return  - compute_MSbar(pars);          
}

int main(){
	
	double pi = 3.1415926535;
	SMDR_Q_in = 173.22; // mass of top-quark, we work here
	float n_of_sigmas = 1;
	
	ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
	
	minimizer->SetMaxFunctionCalls(10000);
	minimizer->SetMaxIterations(1000);
	minimizer->SetTolerance(0.01);
	minimizer->SetPrintLevel(0);

	ROOT::Math::Functor functor(&compute_MSbar, 19);

	minimizer->SetFunction(functor);

	minimizer->SetFixedVariable(0, "Q", 173.22);
	
	// minimizer->SetLimitedVariable(no, "name", initial, step, lower_bound, higher_bound)
	minimizer->SetLimitedVariable(1, "alphaS_5_MZ", SMDR_alphaS_MZ_EXPT, SMDR_alphaS_MZ_EXPT_UNC/10, SMDR_alphaS_MZ_EXPT - n_of_sigmas*SMDR_alphaS_MZ_EXPT_UNC, SMDR_alphaS_MZ_EXPT + n_of_sigmas*SMDR_alphaS_MZ_EXPT_UNC);
	minimizer->SetLimitedVariable(2, "alpha", SMDR_alpha_EXPT, SMDR_alpha_EXPT_UNC/10, SMDR_alpha_EXPT - n_of_sigmas*SMDR_alpha_EXPT_UNC, SMDR_alpha_EXPT + n_of_sigmas*SMDR_alpha_EXPT_UNC);
	minimizer->SetLimitedVariable(3, "GFermi", SMDR_GFermi_EXPT, SMDR_GFermi_EXPT_UNC/10, SMDR_GFermi_EXPT - n_of_sigmas*SMDR_GFermi_EXPT_UNC, SMDR_GFermi_EXPT + n_of_sigmas*SMDR_GFermi_EXPT_UNC);
	minimizer->SetLimitedVariable(4, "MZ", SMDR_MZ_EXPT, SMDR_MZ_EXPT_UNC/10, SMDR_MZ_EXPT - n_of_sigmas*SMDR_MZ_EXPT_UNC, SMDR_MZ_EXPT + n_of_sigmas*SMDR_MZ_EXPT_UNC);
	minimizer->SetLimitedVariable(5, "Mh", SMDR_Mh_EXPT, SMDR_Mh_EXPT_UNC/10, SMDR_Mh_EXPT - n_of_sigmas*SMDR_Mh_EXPT_UNC, SMDR_Mh_EXPT + n_of_sigmas*SMDR_Mh_EXPT_UNC);
	minimizer->SetLimitedVariable(6, "Mt", SMDR_Mt_EXPT, SMDR_Mt_EXPT_UNC/10, SMDR_Mt_EXPT - n_of_sigmas*SMDR_Mt_EXPT_UNC, SMDR_Mt_EXPT + n_of_sigmas*SMDR_Mt_EXPT_UNC);
	minimizer->SetLimitedVariable(7, "mbmb", SMDR_mbmb_EXPT, SMDR_mbmb_EXPT_UNC_hi/5, SMDR_mbmb_EXPT - n_of_sigmas*SMDR_mbmb_EXPT_UNC_lo, SMDR_mbmb_EXPT + n_of_sigmas*SMDR_mbmb_EXPT_UNC_hi);
	minimizer->SetLimitedVariable(8, "Dah5MZ", SMDR_Delta_alpha_had_5_MZ_EXPT, SMDR_Delta_alpha_had_5_MZ_EXPT_UNC/10, SMDR_Delta_alpha_had_5_MZ_EXPT - n_of_sigmas*SMDR_Delta_alpha_had_5_MZ_EXPT_UNC, SMDR_Delta_alpha_had_5_MZ_EXPT + n_of_sigmas*SMDR_Delta_alpha_had_5_MZ_EXPT_UNC);
	
	minimizer->SetFixedVariable(9, "error_tolerance", 1.e-12);
	
	minimizer->SetFixedVariable(10, "loop_config[0]", 1);
	minimizer->SetFixedVariable(11, "loop_config[1]", 1);
	minimizer->SetFixedVariable(12, "loop_config[2]", 1);
	minimizer->SetFixedVariable(13, "loop_config[3]", 1);
	minimizer->SetFixedVariable(14, "loop_config[4]", 1);
	minimizer->SetFixedVariable(15, "loop_config[5]", 1);
	minimizer->SetFixedVariable(16, "loop_config[6]", 1);
	minimizer->SetFixedVariable(17, "loop_config[7]", 1);
	
	minimizer->SetFixedVariable(18, "what", 1);
	//extracting lower bounds
	minimizer->Minimize();
	Double_t g_low = minimizer->MinValue();
	
	minimizer->SetVariableValue(18, 2);
	minimizer->Minimize();
	Double_t gp_low = minimizer->MinValue();
	
	minimizer->SetVariableValue(18, 3);
	minimizer->Minimize();
	Double_t g3_low = minimizer->MinValue();
	
	minimizer->SetVariableValue(18, 4);
	minimizer->Minimize();
	Double_t yt_low = minimizer->MinValue();
	
	minimizer->SetVariableValue(18, 5);
	minimizer->Minimize();
	Double_t yb_low = minimizer->MinValue();
	
	minimizer->SetVariableValue(18, 6);
	minimizer->Minimize();
	Double_t lambda_low = minimizer->MinValue();
	
	// extracting upper bounds
	ROOT::Math::Functor functor2(&compute_MSbar_butminus, 19);

	minimizer->SetFunction(functor2);
	
	minimizer->SetVariableValue(18, 1);
	minimizer->Minimize();
	Double_t g_high = - minimizer->MinValue();
	
	minimizer->SetVariableValue(18, 2);
	minimizer->Minimize();
	Double_t gp_high = - minimizer->MinValue();
	
	minimizer->SetVariableValue(18, 3);
	minimizer->Minimize();
	Double_t g3_high = - minimizer->MinValue();
	
	minimizer->SetVariableValue(18, 4);
	minimizer->Minimize();
	Double_t yt_high = - minimizer->MinValue();
	
	minimizer->SetVariableValue(18, 5);
	minimizer->Minimize();
	Double_t yb_high = - minimizer->MinValue();
	
	minimizer->SetVariableValue(18, 6);
	minimizer->Minimize();
	Double_t lambda_high = - minimizer->MinValue();
	
	cout << endl;
	
	cout << "g is (" << g_low << " ... " << g_high <<")" << endl;
	cout << "g' is (" << gp_low << " ... " << gp_high <<")" << endl;
	cout << "g3 is (" << g3_low << " ... " << g3_high <<")" << endl;
	cout << "yt is (" << yt_low << " ... " << yt_high <<")" << endl;
	cout << "yb is (" << yb_low << " ... " << yb_high <<")" << endl;
	cout << "lambda is (" << lambda_low << " ... " << lambda_high <<")" << endl;
	
	cout << endl;
	
	cout << "x1 is (" << g_low*g_low*5 / (48*pi*pi) << " ... " << g_high*g_high*5 / (48*pi*pi) <<")" << endl;
	cout << "x2 is (" << gp_low*gp_low / (16*pi*pi) << " ... " << gp_high*gp_high / (16*pi*pi) <<")" << endl;
	cout << "x3 is (" << g3_low*g3_low / (16*pi*pi) << " ... " << g3_high*g3_high / (16*pi*pi) <<")" << endl;
	cout << "y1 is (" << yt_low*yt_low / (16*pi*pi) << " ... " << yt_high*yt_high / (16*pi*pi) <<")" << endl;
	cout << "y2 is (" << yb_low*yb_low / (16*pi*pi) << " ... " << yb_high*yb_high / (16*pi*pi) <<")" << endl;
	cout << "z is (" << lambda_low / (16*pi*pi) << " ... " << lambda_high / (16*pi*pi) <<")" << endl;

	
	return 0;
}
	
