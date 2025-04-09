// old program that does not use Fit_Inputs(), but rather uses Minuit to fit. makes incorrect results and very slow, so the idea was abandoned
// to launch the program: 
// g++ `root-config --cflags` direct_minimization.cpp `root-config --libs` -lm -lsmdr -ltsil -l3vil

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

// experimental data written explicitly so it could be handled in this script
SMDR_REAL SMDR_alphaS_5_MZ_Exp = 0.1179; 
SMDR_REAL SMDR_alphaS_5_MZ_Exp_UNC = 0.0010;// SMDR suggests this??..
SMDR_REAL SMDR_alpha_Exp = 0.0072973525693;
SMDR_REAL SMDR_alpha_Exp_UNC = 0.0000000000011;
SMDR_REAL SMDR_GFermi_Exp = 0.000011663787;
SMDR_REAL SMDR_GFermi_Exp_UNC = 0.000000000006;
SMDR_REAL SMDR_MZ_Exp = 91.1876; 
SMDR_REAL SMDR_MZ_Exp_UNC = 0.0021;
SMDR_REAL SMDR_Mh_Exp = 125.25;
SMDR_REAL SMDR_Mh_Exp_UNC = 0.17;
SMDR_REAL SMDR_Mt_Exp = 172.5;
SMDR_REAL SMDR_Mt_Exp_UNC = 0.7;
SMDR_REAL SMDR_mbmb_Exp = 4.18; 
SMDR_REAL SMDR_mbmb_Exp_UNC = 0.03;
SMDR_REAL SMDR_Delta_alpha_had_5_MZ_Exp = 0.02766;
SMDR_REAL SMDR_Delta_alpha_had_5_MZ_Exp_UNC = 0.00007; // from smdr_pdg_2021.h, there from Review of Particle Properties 2020 

Double_t RSS(const Double_t* pars){ // sum of squares
	// pars of the form {gp,g,g3,yt,yb,lambda,m2,Lambda,Deltaalpha, loop_config[0]...}
	
	SMDR_Q_in = 173.22; // mass of top-quark, we work here
	
	SMDR_lambda_in = pars[5]; 
	SMDR_m2_in = pars[6]; 
	SMDR_g3_in = pars[2];
	SMDR_g_in = pars[1];
	SMDR_gp_in = pars[0];
	SMDR_yt_in = pars[3];
	SMDR_yb_in = pars[4]; 
	SMDR_yc_in = 0;
	SMDR_ys_in = 0;
	SMDR_yu_in = 0;
	SMDR_yd_in = 0;
	SMDR_ytau_in = 0;
	SMDR_ymu_in = 0;
	SMDR_ye_in = 0; // disabled alltogether
	SMDR_Delta_alpha_had_5_MZ_in = pars[7];
	SMDR_v_in = SMDR_Eval_vev (173.22, 1);
	
	SMDR_Load_Inputs ();
	
	SMDR_Eval_Mt_pole (173.22, 1, pars[8], pars[9], &SMDR_Mt_pole, &SMDR_Gammat_pole);
	SMDR_Eval_Mh_pole (173.22, pars[10], &SMDR_Mh_pole, &SMDR_Gammah_pole);
	SMDR_Eval_MZ_pole (173.22, pars[11], &SMDR_MZ_pole, &SMDR_GammaZ_pole,
		       &SMDR_MZ_PDG, &SMDR_GammaZ_PDG);

	/* This one is computed only because it is needed by SMDR_Eval_Gauge() */
	SMDR_Eval_MW_pole (173.22, pars[12], &SMDR_MW_pole, &SMDR_GammaW_pole,
		       &SMDR_MW_PDG, &SMDR_GammaW_PDG);

	SMDR_GFermi = SMDR_Eval_GFermi (SMDR_Mt_pole, pars[13]);
	SMDR_Eval_Gauge (SMDR_Mt_pole, SMDR_Mh_pole, SMDR_MW_pole);
	SMDR_Eval_QCDQED_at_MZ (SMDR_MZ_PDG, 173.22, pars[14]);
	// here i should use Q or MZ???????
	SMDR_mbmb = SMDR_Eval_mbmb(173.22, pars[15]);
	
	SMDR_Load_Inputs ();
	
	return (SMDR_alphaS_5_MZ - SMDR_alphaS_5_MZ_Exp)*(SMDR_alphaS_5_MZ - SMDR_alphaS_5_MZ_Exp)/(SMDR_alphaS_5_MZ_Exp_UNC*SMDR_alphaS_5_MZ_Exp_UNC) + (SMDR_alpha - SMDR_alpha_Exp)*(SMDR_alpha - SMDR_alpha_Exp)/(SMDR_alpha_Exp_UNC*SMDR_alpha_Exp_UNC) + (SMDR_GFermi - SMDR_GFermi_Exp)*(SMDR_GFermi - SMDR_GFermi_Exp)/(SMDR_GFermi_Exp_UNC*SMDR_GFermi_Exp_UNC) + (SMDR_MZ_PDG - SMDR_MZ_Exp)*(SMDR_MZ_PDG - SMDR_MZ_Exp)/(SMDR_MZ_Exp_UNC*SMDR_MZ_Exp_UNC) + (SMDR_Mh_pole - SMDR_Mh_Exp)*(SMDR_Mh_pole - SMDR_Mh_Exp)/(SMDR_Mh_Exp_UNC*SMDR_Mh_Exp_UNC) + (SMDR_Mt_pole - SMDR_Mt_Exp)*(SMDR_Mt_pole - SMDR_Mt_Exp)/(SMDR_Mt_Exp_UNC*SMDR_Mt_Exp_UNC) + (SMDR_mbmb - SMDR_mbmb_Exp)*(SMDR_mbmb - SMDR_mbmb_Exp)/(SMDR_mbmb_Exp_UNC*SMDR_mbmb_Exp_UNC) + (SMDR_Delta_alpha_had_5_MZ - SMDR_Delta_alpha_had_5_MZ_Exp)*(SMDR_Delta_alpha_had_5_MZ - SMDR_Delta_alpha_had_5_MZ_Exp)/(SMDR_Delta_alpha_had_5_MZ_Exp_UNC*SMDR_Delta_alpha_had_5_MZ_Exp_UNC);

}

int main(){

	SMDR_Q_in = 173.22; // mass of top-quark, we work here
	
	ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit", "Minos");
	
	minimizer->SetMaxFunctionCalls(1000000);
	minimizer->SetMaxIterations(100000);
	minimizer->SetTolerance(0.0001);
	minimizer->SetPrintLevel(1);

	ROOT::Math::Functor functor(&RSS, 16);

	minimizer->SetFunction(functor);

	minimizer->SetLimitedVariable(0, "gp", 0.35838, 0.00001, 0.33, 0.37);
	minimizer->SetLimitedVariable(1, "g", 0.64812, 0.00001, 0.62, 0.66);
	minimizer->SetLimitedVariable(2, "g3", 1.1654, 0.0001, 1.15, 1.17);
	minimizer->SetLimitedVariable(3, "yt", 0.93517, 0.00001, 0.92, 0.95);
	minimizer->SetLimitedVariable(4, "yb", 0.01706, 0.0001, 0.0155, 0.0185);
	minimizer->SetLimitedVariable(5, "lambda", 0.1274, 0.0001, 0.126, 0.129);
	minimizer->SetLimitedVariable(6, "m2", -8624, 10, -8400, -8900);
	//minimizer->SetLimitedVariable(7, "v", 246.2, 0.1, 230, 260);
	minimizer->SetLimitedVariable(7, "Dah5MZ", 0.02766, 0.0001, 0.025, 0.029); // Delta_alpha_had_5_MZ
	
	minimizer->SetFixedVariable(8, "loop_config[0]", 1);
	minimizer->SetFixedVariable(9, "loop_config[1]", 1);
	minimizer->SetFixedVariable(10, "loop_config[2]", 1);
	minimizer->SetFixedVariable(11, "loop_config[3]", 1);
	minimizer->SetFixedVariable(12, "loop_config[4]", 1);
	minimizer->SetFixedVariable(13, "loop_config[5]", 1);
	minimizer->SetFixedVariable(14, "loop_config[6]", 1);
	minimizer->SetFixedVariable(15, "loop_config[7]", 1);
	

	minimizer->Minimize();

	
	return 0;
}
