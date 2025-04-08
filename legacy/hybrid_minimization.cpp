// first program to somewhat successfully extract parametric errors by minimizing the corresponding MSbar within the "cube of errors" of on-shell's. actually, gives pretty close results, but anyway is NOT rigorous
// to launch the program: 
// g++ `root-config --cflags` hybrid_minimization.cpp `root-config --libs` -lm -lsmdr -ltsil -l3vil

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

#include "../loop_Fit_Inputs.cpp"
#include "../smdr_pdg_2025.h"

float n_of_sigmas = 1;

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

void init_minimizer(ROOT::Math::Minimizer* minimizer)
{
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
	
	minimizer->SetFixedVariable(9, "error_tolerance", 1.e-6);
	
	minimizer->SetFixedVariable(10, "loop_config[0]", 1);
	minimizer->SetFixedVariable(11, "loop_config[1]", 1);
	minimizer->SetFixedVariable(12, "loop_config[2]", 1);
	minimizer->SetFixedVariable(13, "loop_config[3]", 1);
	minimizer->SetFixedVariable(14, "loop_config[4]", 1);
	minimizer->SetFixedVariable(15, "loop_config[5]", 1);
	minimizer->SetFixedVariable(16, "loop_config[6]", 1);
	minimizer->SetFixedVariable(17, "loop_config[7]", 1);
}

int main(){
	
	double pi = 3.1415926535;
	SMDR_Q_in = 173.22; // mass of top-quark, we work here
	
	ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
	
	minimizer->SetMaxFunctionCalls(10000);
	minimizer->SetMaxIterations(10000);
	minimizer->SetTolerance(0.001);
	minimizer->SetPrintLevel(0);

	ROOT::Math::Functor functor(&compute_MSbar, 19);

	minimizer->SetFunction(functor);

	init_minimizer(minimizer);
	
	minimizer->SetFixedVariable(18, "what", 1);
	//extracting lower bounds
	minimizer->Minimize();
	Double_t g_low = minimizer->MinValue();
	
	minimizer->SetVariableValue(18, 2);
	init_minimizer(minimizer);
	minimizer->Minimize();
	Double_t gp_low = minimizer->MinValue();
	
	minimizer->SetVariableValue(18, 3);
	init_minimizer(minimizer);
	minimizer->Minimize();
	Double_t g3_low = minimizer->MinValue();
	
	minimizer->SetVariableValue(18, 4);
	init_minimizer(minimizer);
	minimizer->Minimize();
	Double_t yt_low = minimizer->MinValue();
	
	minimizer->SetVariableValue(18, 5);
	init_minimizer(minimizer);
	minimizer->Minimize();
	Double_t yb_low = minimizer->MinValue();
	
	minimizer->SetVariableValue(18, 6);
	init_minimizer(minimizer);
	minimizer->Minimize();
	Double_t lambda_low = minimizer->MinValue();
	
	// extracting upper bounds
	ROOT::Math::Functor functor2(&compute_MSbar_butminus, 19);

	minimizer->SetFunction(functor2);
	
	init_minimizer(minimizer);
	
	minimizer->SetVariableValue(18, 1);
	minimizer->Minimize();
	Double_t g_high = - minimizer->MinValue();
	
	minimizer->SetVariableValue(18, 2);
	init_minimizer(minimizer);
	minimizer->Minimize();
	Double_t gp_high = - minimizer->MinValue();
	
	minimizer->SetVariableValue(18, 3);
	init_minimizer(minimizer);
	minimizer->Minimize();
	Double_t g3_high = - minimizer->MinValue();
	
	minimizer->SetVariableValue(18, 4);
	init_minimizer(minimizer);
	minimizer->Minimize();
	Double_t yt_high = - minimizer->MinValue();
	
	minimizer->SetVariableValue(18, 5);
	init_minimizer(minimizer);
	minimizer->Minimize();
	Double_t yb_high = - minimizer->MinValue();
	
	minimizer->SetVariableValue(18, 6);
	init_minimizer(minimizer);
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
	
