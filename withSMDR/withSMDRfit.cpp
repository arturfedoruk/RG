// to launch the program: 
// g++ `root-config --cflags` withSMDRfit.cpp `root-config --libs` -lm -lsmdr -ltsil -l3vil

#include "smdr.h"
#include "iostream"
#include "fstream"
#include "string"
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TError.h"
using namespace std;

TFile *file;
TTree *tree;
TTree *tree_assess;

string filename_default = "net_data.root"; // net_datas/net_333221.root
string treename_default = "tree_222222_5,_:5_5_5_5_5_5_5"; // tree_333221_5,_:5_5_5_5_5_5_5
string treename_assess_default = "tree_222222_5,_:5_5_5_5_5_5_5";

double alphaS, Mt, Mh, MZ, Delta, mbmb, GF, g, gp, g3, yt, yb, lambda;

// estimated sigmas in order to pass them into chi^2 in minimization (estimated based on 222222 random generation)
double lambda_sigma = 0.00023;
double g_sigma = 2.5e-5;
double gp_sigma = 1.7e-5;
double g3_sigma = 0.004;
double yt_sigma = 0.0017;
double yb_sigma = 6.7e-5;

Double_t fit_lambda(Double_t* xx, const Double_t* par){
	double deltah = (xx[0] - SMDR_Mh_EXPT)/SMDR_Mh_EXPT;
	double deltat = (xx[1] - SMDR_Mt_EXPT)/SMDR_Mt_EXPT;
	double deltaZ = (xx[2] - SMDR_MZ_EXPT)/SMDR_MZ_EXPT;
	double deltaS = (xx[3] - SMDR_alphaS_MZ_EXPT)/SMDR_alphaS_MZ_EXPT;
	double deltaa = (xx[4] - SMDR_Delta_alpha_had_5_MZ_EXPT)/SMDR_Delta_alpha_had_5_MZ_EXPT;
	double deltab = (xx[5] - SMDR_mbmb_EXPT)/SMDR_mbmb_EXPT;
	double deltaGF = (xx[6] - SMDR_GFermi_EXPT)/ SMDR_GFermi_EXPT;
	return par[15] * (1 + par[0]*deltah + par[1]*deltat + par[2]*deltaZ + par[3]*deltaS + par[4]*deltat*deltat + par[5]*deltat*deltaS + par[6]*deltah*deltah + par[7]*deltah*deltat + par[8]*deltaS*deltaS + par[9]*deltah*deltaS + par[10]*deltat*deltat*deltat + par[11]*deltat*deltat*deltaS + par[12]*deltaa + par[13]*deltab + par[14]*deltaGF);
}

Double_t sum_squares_lambda(const Double_t* pars){
	Double_t sum = 0;
	Double_t onshell[7];
	for(int ientry = 0; ientry < tree->GetEntries(); ientry++){
		tree->GetEntry(ientry);
		onshell[0] = Mh;
		onshell[1] = Mt;
		onshell[2] = MZ;
		onshell[3] = alphaS;
		onshell[4] = Delta;
		onshell[5] = mbmb;
		onshell[6] = GF;
		sum += ( lambda - fit_lambda(onshell, pars) )*( lambda - fit_lambda(onshell, pars) ) / (lambda_sigma*lambda_sigma);
	}
	return sum;
}

Double_t fit_g(Double_t* xx, const Double_t* par){
	double deltah = (xx[0] - SMDR_Mh_EXPT)/SMDR_Mh_EXPT;
	double deltat = (xx[1] - SMDR_Mt_EXPT)/SMDR_Mt_EXPT;
	double deltaZ = (xx[2] - SMDR_MZ_EXPT)/SMDR_MZ_EXPT;
	double deltaS = (xx[3] - SMDR_alphaS_MZ_EXPT)/SMDR_alphaS_MZ_EXPT;
	double deltaa = (xx[4] - SMDR_Delta_alpha_had_5_MZ_EXPT)/SMDR_Delta_alpha_had_5_MZ_EXPT;
	double deltaGF = (xx[6] - SMDR_GFermi_EXPT)/ SMDR_GFermi_EXPT;
	return par[11] * (1 + par[0]*deltah + par[1]*deltat + par[2]*deltaZ + par[3]*deltaS + par[4]*deltat*deltat + par[5]*deltat*deltaS + par[6]*deltaS*deltaS + par[7]*deltaZ*deltaZ + par[8]*deltaZ*deltat + par[9]*deltaa + par[10]*deltaGF);
}

Double_t sum_squares_g(const Double_t* pars){
	Double_t sum = 0;
	Double_t onshell[7];
	for(int ientry = 0; ientry < tree->GetEntries(); ientry++){
		tree->GetEntry(ientry);		
		onshell[0] = Mh;
		onshell[1] = Mt;
		onshell[2] = MZ;
		onshell[3] = alphaS;
		onshell[4] = Delta;
		onshell[5] = mbmb;
		onshell[6] = GF;
		sum += ( g - fit_g(onshell, pars) )*( g - fit_g(onshell, pars) ) /(g_sigma*g_sigma);	
	}	
	return sum;
}

Double_t fit_gp(Double_t* xx, const Double_t* par){
	double deltah = (xx[0] - SMDR_Mh_EXPT)/SMDR_Mh_EXPT;
	double deltat = (xx[1] - SMDR_Mt_EXPT)/SMDR_Mt_EXPT;
	double deltaZ = (xx[2] - SMDR_MZ_EXPT)/SMDR_MZ_EXPT;
	double deltaS = (xx[3] - SMDR_alphaS_MZ_EXPT)/SMDR_alphaS_MZ_EXPT;
	double deltaa = (xx[4] - SMDR_Delta_alpha_had_5_MZ_EXPT)/SMDR_Delta_alpha_had_5_MZ_EXPT;
	double deltaGF = (xx[6] - SMDR_GFermi_EXPT)/ SMDR_GFermi_EXPT;
	return par[9] * (1 + par[0]*deltah + par[1]*deltat + par[2]*deltaZ + par[3]*deltaS + par[4]*deltat*deltat + par[5]*deltaZ*deltat + par[6]*deltaZ*deltaZ + par[7]*deltaa + par[8]*deltaGF);
}

Double_t sum_squares_gp(const Double_t* pars){
	Double_t sum = 0;
	Double_t onshell[7];	
	for(int ientry = 0; ientry < tree->GetEntries(); ientry++){		
		tree->GetEntry(ientry);		
		onshell[0] = Mh;
		onshell[1] = Mt;
		onshell[2] = MZ;
		onshell[3] = alphaS;
		onshell[4] = Delta;
		onshell[5] = mbmb;
		onshell[6] = GF;
		sum += ( gp - fit_gp(onshell, pars) )*( gp - fit_gp(onshell, pars) ) / (gp_sigma*gp_sigma);
	}	
	return sum;
}

Double_t fit_g3(Double_t* xx, const Double_t* par){
	double deltah = (xx[0] - SMDR_Mh_EXPT)/SMDR_Mh_EXPT;
	double deltat = (xx[1] - SMDR_Mt_EXPT)/SMDR_Mt_EXPT;
	double deltaZ = (xx[2] - SMDR_MZ_EXPT)/SMDR_MZ_EXPT;
	double deltaS = (xx[3] - SMDR_alphaS_MZ_EXPT)/SMDR_alphaS_MZ_EXPT;
	double deltaa = (xx[4] - SMDR_Delta_alpha_had_5_MZ_EXPT)/SMDR_Delta_alpha_had_5_MZ_EXPT;
	double deltab = (xx[5] - SMDR_mbmb_EXPT)/SMDR_mbmb_EXPT;
	double deltaGF = (xx[6] - SMDR_GFermi_EXPT)/ SMDR_GFermi_EXPT;
	return par[8] * (1 + par[0]*deltah + par[1]*deltat + par[2]*deltaZ + par[3]*deltaS + par[4]*deltat*deltat + par[5]*deltat*deltaS + par[6]*deltaS*deltaS + par[7]*deltaS*deltaS*deltaS);
}

Double_t sum_squares_g3(const Double_t* pars){
	Double_t sum = 0;
	Double_t onshell[7];	
	for(int ientry = 0; ientry < tree->GetEntries(); ientry++){		
		tree->GetEntry(ientry);
		onshell[0] = Mh;
		onshell[1] = Mt;
		onshell[2] = MZ;
		onshell[3] = alphaS;
		onshell[4] = Delta;
		onshell[5] = mbmb;
		onshell[6] = GF;
		sum += ( g3 - fit_g3(onshell, pars) )*( g3 - fit_g3(onshell, pars) ) / (g3_sigma*g3_sigma);	
	}	
	return sum;
}

Double_t fit_yt(Double_t* xx, const Double_t* par){
	double deltah = (xx[0] - SMDR_Mh_EXPT)/SMDR_Mh_EXPT;
	double deltat = (xx[1] - SMDR_Mt_EXPT)/SMDR_Mt_EXPT;
	double deltaZ = (xx[2] - SMDR_MZ_EXPT)/SMDR_MZ_EXPT;
	double deltaS = (xx[3] - SMDR_alphaS_MZ_EXPT)/SMDR_alphaS_MZ_EXPT;
	double deltaa = (xx[4] - SMDR_Delta_alpha_had_5_MZ_EXPT)/SMDR_Delta_alpha_had_5_MZ_EXPT;
	return par[9] * (1 + par[0]*deltah + par[1]*deltat + par[2]*deltaZ + par[3]*deltaS + par[4]*deltat*deltat + par[5]*deltat*deltaS + par[6]*deltaS*deltaS + par[7]*deltah*deltat + par[8]*deltaa);
}

Double_t sum_squares_yt(const Double_t* pars){
	Double_t sum = 0;
	Double_t onshell[7];	
	for(int ientry = 0; ientry < tree->GetEntries(); ientry++){		
		tree->GetEntry(ientry);		
		onshell[0] = Mh;
		onshell[1] = Mt;
		onshell[2] = MZ;
		onshell[3] = alphaS;
		onshell[4] = Delta;
		onshell[5] = mbmb;
		onshell[6] = GF;
		sum += ( yt - fit_yt(onshell, pars) )*( yt - fit_yt(onshell, pars) ) / (yt_sigma*yt_sigma);
	}	
	return sum;
}


Double_t fit_yb(Double_t* xx, const Double_t* par){
	double deltah = (xx[0] - SMDR_Mh_EXPT)/SMDR_Mh_EXPT;
	double deltat = (xx[1] - SMDR_Mt_EXPT)/SMDR_Mt_EXPT;
	double deltaZ = (xx[2] - SMDR_MZ_EXPT)/SMDR_MZ_EXPT;
	double deltaS = (xx[3] - SMDR_alphaS_MZ_EXPT)/SMDR_alphaS_MZ_EXPT;
	double deltab = (xx[5] - SMDR_mbmb_EXPT)/SMDR_mbmb_EXPT;
	return par[10] * (1 + par[0]*deltab + par[1]*deltah + par[2]*deltat + par[3]*deltaZ + par[4]*deltaS + par[5]*deltab*deltab + par[6]*deltab*deltaS + par[7]*deltaS*deltaS + par[8]*deltat*deltaS + par[9]*deltaS*deltaS*deltaS);
}

Double_t sum_squares_yb(const Double_t* pars){
	Double_t sum = 0;
	Double_t onshell[7];	
	for(int ientry = 0; ientry < tree->GetEntries(); ientry++){		
		tree->GetEntry(ientry);
		onshell[0] = Mh;
		onshell[1] = Mt;
		onshell[2] = MZ;
		onshell[3] = alphaS;
		onshell[4] = Delta;
		onshell[5] = mbmb;
		onshell[6] = GF;
		sum += ( yb - fit_yb(onshell, pars) )*( yb - fit_yb(onshell, pars) ) / (yb_sigma*yb_sigma);
	}
	return sum;
}

//------------------------------------------------------------------------------
//-----------------MAIN---------------------------------------------------------
int main(){

	#include "smdr_pdg_2025.h"	

	char if_default_files;
	cout << "Use default files? y/n" << endl;
	cin >> if_default_files;
	string filename_temp, treename_temp, treename_assess_temp;
	if (if_default_files == 'n'){
		cout << "Enter file name " << endl;
		cin >> filename_temp;
		cout << "Enter name of tree to learn " << endl;
		cin >> treename_temp;
		cout << "Enter name of tree to assess " << endl;
		cin >> treename_assess_temp;
	} else {
		filename_temp = filename_default;
		treename_temp = treename_default;
		treename_assess_temp = treename_assess_default;
	}
	
	ofstream outfile(("fit_outputs/output_" + treename_temp + ".txt").c_str(), ios::app | ios::out);
	
	const char* filename = filename_temp.c_str();
	const char* treename = treename_temp.c_str();
	const char* treename_assess = treename_assess_temp.c_str();
	
	file = new TFile(filename);
	tree = (TTree*)file->Get(treename);
	tree_assess = (TTree*)file->Get(treename_assess);
	
	if (!tree) cout << "LEARNING TREE NOT FOUND" << endl;
	if (!tree_assess) cout << "ASSESS TREE NOT FOUND" << endl;
	
	// Fitting
	
	cout << endl << "FITTING\n" << endl;	
	
	tree->SetBranchAddress("alphaS_MZ",&alphaS);
	tree->SetBranchAddress("Mt",&Mt);
	tree->SetBranchAddress("Mh",&Mh);
	tree->SetBranchAddress("MZ",&MZ);
	tree->SetBranchAddress("Delta_alpha",&Delta);
	tree->SetBranchAddress("mbmb",&mbmb);
	tree->SetBranchAddress("GFermi",&GF);
		
	tree_assess->SetBranchAddress("alphaS_MZ",&alphaS);
	tree_assess->SetBranchAddress("Mt",&Mt);
	tree_assess->SetBranchAddress("Mh",&Mh);
	tree_assess->SetBranchAddress("MZ",&MZ);
	tree_assess->SetBranchAddress("Delta_alpha",&Delta);
	tree_assess->SetBranchAddress("mbmb",&mbmb);
	tree_assess->SetBranchAddress("GFermi",&GF);
	
	int what;
	cout << "What are we fitting tonight? \n 1 - g, 2 - g', 3 - g3, 4 - yt, 5 - yb, 6 - lambda" << endl;
	cin >> what;
	
	if(what==6)
	{	
		tree->SetBranchAddress("lambda",&lambda);

		ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");			
		minimizer->SetMaxFunctionCalls(10000);
		minimizer->SetMaxIterations(10000);
		minimizer->SetTolerance(0.000001);
		minimizer->SetPrintLevel(1);

		ROOT::Math::Functor functor(&sum_squares_lambda, 16); // <---!
		minimizer->SetFunction(functor);

		minimizer->SetVariable(0, "c h", 1.68e-3, 1e-5);
		minimizer->SetVariable(1, "c t", -1.49e-4, 1e-6);
		minimizer->SetVariable(2, "c Z", -3.5e-7, 1e-8);
		minimizer->SetVariable(3, "c S", -2.2e-7, 1e-8);
		minimizer->SetVariable(4, "c tt", 1.5e-5, 1e-6);
		minimizer->SetVariable(5, "c tS", 4e-6, 1e-7);
		minimizer->SetVariable(6, "c hh", 7e-7, 1e-8);
		minimizer->SetVariable(7, "c ht", -6.1e-7, 1e-8);
		minimizer->SetVariable(8, "c SS", 3e-7, 1e-8);
		minimizer->SetVariable(9, "c hS", 6.4e-8, 1e-9);
		minimizer->SetVariable(10, "c ttt", 1.9e-7, 1e-8);
		minimizer->SetVariable(11, "c ttS", -7.6e-8, 1e-9);
		minimizer->SetVariable(12, "c a", 3.4e-7, 1e-8);			
		minimizer->SetVariable(13, "c b", 4.5e-5, 1e-6);
		minimizer->SetVariable(14, "c GF", 0.95, 0.01);			
		minimizer->SetVariable(15, "lambda 0", 0.127, 3e-5);
		
		minimizer->Minimize();
		const double* fit_results = minimizer->X();
		
		double lambda_0 = fit_results[15];
		
		// Propogation of errors
		
		double Delta_lambda = lambda_0 * TMath::Sqrt( fit_results[0]*fit_results[0]*SMDR_Mh_EXPT_UNC*SMDR_Mh_EXPT_UNC/SMDR_Mh_EXPT/SMDR_Mh_EXPT + fit_results[1]*fit_results[1]*SMDR_Mt_EXPT_UNC*SMDR_Mt_EXPT_UNC/SMDR_Mt_EXPT/SMDR_Mt_EXPT + fit_results[2]*fit_results[2]*SMDR_MZ_EXPT_UNC*SMDR_MZ_EXPT_UNC/SMDR_MZ_EXPT/SMDR_MZ_EXPT + fit_results[3]*fit_results[3]*SMDR_alphaS_MZ_EXPT_UNC*SMDR_alphaS_MZ_EXPT_UNC/SMDR_alphaS_MZ_EXPT/SMDR_alphaS_MZ_EXPT + fit_results[12]*fit_results[12]*SMDR_Delta_alpha_had_5_MZ_EXPT_UNC*SMDR_Delta_alpha_had_5_MZ_EXPT_UNC/SMDR_Delta_alpha_had_5_MZ_EXPT/SMDR_Delta_alpha_had_5_MZ_EXPT + fit_results[13]*fit_results[13]*SMDR_mbmb_EXPT_UNC_hi*SMDR_mbmb_EXPT_UNC_hi/SMDR_mbmb_EXPT/SMDR_mbmb_EXPT + fit_results[14]*fit_results[14]*SMDR_GFermi_EXPT_UNC*SMDR_GFermi_EXPT_UNC/SMDR_GFermi_EXPT/SMDR_GFermi_EXPT );
		//double Delta_lambda = lambda_0 * ( TMath::Abs(fit_results[0])*SMDR_Mh_EXPT_UNC/SMDR_Mh_EXPT + TMath::Abs(fit_results[1])*SMDR_Mt_EXPT_UNC/SMDR_Mt_EXPT + TMath::Abs(fit_results[2])*SMDR_MZ_EXPT_UNC/SMDR_MZ_EXPT + TMath::Abs(fit_results[3])*SMDR_alphaS_MZ_EXPT_UNC/SMDR_alphaS_MZ_EXPT + TMath::Abs(fit_results[12])*SMDR_Delta_alpha_had_5_MZ_EXPT_UNC/SMDR_Delta_alpha_had_5_MZ_EXPT + TMath::Abs(fit_results[13])*SMDR_mbmb_EXPT_UNC_hi/SMDR_mbmb_EXPT + TMath::Abs(fit_results[14])*SMDR_GFermi_EXPT_UNC/SMDR_GFermi_EXPT );
		
		cout << "\n++++++++++++++++++++++++++++++++++++++\nPropogation of errors: \nDelta lambda = " << Delta_lambda << "\nDelta z = " << Delta_lambda / (16. * TMath::Pi()*TMath::Pi()) << endl;
		
		outfile << "\n++++++++++++++++++++++++++++++++++++++\nPropogation of errors: \nDelta lambda = " << Delta_lambda << "\nDelta z = " << Delta_lambda / (16. * TMath::Pi()*TMath::Pi()) << endl;
		
		// Assessing
		
		cout << endl << "\n++++++++++++++++++++++++++++++++++++++\nASSESSING\n" << endl;
		tree_assess->SetBranchAddress("lambda",&lambda);
		
		double lambda_min = lambda_0;
		double lambda_max = lambda_0;
		double lambda_errmean = 0;
		double lambda_errmax = 0;
		
		double lambda_fromfit, lambda_err;
		
		int N = tree_assess->GetEntries();
		
		double onshell[7];
		
		for(int ientry = 0; ientry < N; ientry++){
			
			tree_assess->GetEntry(ientry);
			
			if(lambda < lambda_min) lambda_min = lambda;
			if(lambda > lambda_max) lambda_max = lambda;
			
			onshell[0] = Mh;
			onshell[1] = Mt;
			onshell[2] = MZ;
			onshell[3] = alphaS;
			onshell[4] = Delta;
			onshell[5] = mbmb;
			onshell[6] = GF;
			
			lambda_fromfit = fit_lambda(onshell, fit_results);
			
			lambda_err = TMath::Abs( lambda_fromfit - lambda );
			
			lambda_errmean += lambda_err / N;
			
			if (lambda_err > lambda_errmax) lambda_errmax = lambda_err;
			
			cout << "Processed event: " << ientry+1 << "/" << N << "\r";
			cout.flush();
		
		}
		
		cout << "\n\nlambda_mean = " << lambda_0 << "\nlambda_min = " << lambda_min << "\nlambda_max = " << lambda_max << "\nsigma = " << Delta_lambda << " \nParametric error: " << Delta_lambda / lambda_0 << "\n\nlambda_errmean = " << lambda_errmean << "\nlambda_errmax = " << lambda_errmax << "\nFit mean error: " << lambda_errmean / lambda_0 << "\nFit max error: " << lambda_errmax / lambda_0 << endl; 
		
		outfile << "\nParametric error: " << Delta_lambda / lambda_0 << "\nlambda_errmean = " << lambda_errmean << "\nlambda_errmax = " << lambda_errmax << "\nFit mean error: " << lambda_errmean / lambda_0 << "\nFit max error: " << lambda_errmax / lambda_0 << endl; 
		
		cout << "\n++++++++++++++++++++++++++++++++++++++\nFor table: \n{ " << fit_results[15];
		for (int k = 0; k < 15; k++) cout << ", " << fit_results[k];
		cout << ", " << Delta_lambda / lambda_0 << ", " << lambda_errmean / lambda_0 << ", " << lambda_errmax / lambda_0 << " }" << endl;
		
		outfile << "\nFor table: \n{ " << fit_results[15];
		for (int k = 0; k < 15; k++) outfile << ", " << fit_results[k];
		outfile << ", " << Delta_lambda / lambda_0 << ", " << lambda_errmean / lambda_0 << ", " << lambda_errmax / lambda_0 << " }\n++++++++++++++++++++++++++++++++++++++" << endl;
		
	}
	
	if(what==1)
	{
		
		tree->SetBranchAddress("g",&g);

		ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");			
		minimizer->SetMaxFunctionCalls(10000);
		minimizer->SetMaxIterations(10000);
		minimizer->SetTolerance(0.000001);
		minimizer->SetPrintLevel(1);

		ROOT::Math::Functor functor(&sum_squares_g, 12); // <---!
		minimizer->SetFunction(functor);

		minimizer->SetVariable(0, "c h", -8.5e-7, 1e-8);
		minimizer->SetVariable(1, "c t", 5.735e-5, 1e-8);
		minimizer->SetVariable(2, "c Z", 1.558e-5, 1e-8);
		minimizer->SetVariable(3, "c S", -5.97e-6, 1e-8);
		minimizer->SetVariable(4, "c tt", 1.9e-7, 1e-8);
		minimizer->SetVariable(5, "c tS", -7.8e-8, 1e-9);
		minimizer->SetVariable(6, "c SS", 5e-9, 1e-10);
		minimizer->SetVariable(7, "c ZZ", -1e-10, 1e-11);
		minimizer->SetVariable(8, "c Zt", -5e-10, 1e-11);
		minimizer->SetVariable(9, "c a", -2.3e-5, 1e-6);
		minimizer->SetVariable(10, "c GF", 0.71, 0.01);			
		minimizer->SetVariable(11, "g 0", 0.6474, 1e-5);
		
		minimizer->Minimize();
		const double* fit_results = minimizer->X(); 
		
		double g_0 = fit_results[11];
		
		//Propogation of errors
		
		double Delta_g = g_0 * TMath::Sqrt( fit_results[0]*fit_results[0]*SMDR_Mh_EXPT_UNC*SMDR_Mh_EXPT_UNC/SMDR_Mh_EXPT/SMDR_Mh_EXPT + fit_results[1]*fit_results[1]*SMDR_Mt_EXPT_UNC*SMDR_Mt_EXPT_UNC/SMDR_Mt_EXPT/SMDR_Mt_EXPT + fit_results[2]*fit_results[2]*SMDR_MZ_EXPT_UNC*SMDR_MZ_EXPT_UNC/SMDR_MZ_EXPT/SMDR_MZ_EXPT + fit_results[3]*fit_results[3]*SMDR_alphaS_MZ_EXPT_UNC*SMDR_alphaS_MZ_EXPT_UNC/SMDR_alphaS_MZ_EXPT/SMDR_alphaS_MZ_EXPT + fit_results[9]*fit_results[9]*SMDR_Delta_alpha_had_5_MZ_EXPT_UNC*SMDR_Delta_alpha_had_5_MZ_EXPT_UNC/SMDR_Delta_alpha_had_5_MZ_EXPT/SMDR_Delta_alpha_had_5_MZ_EXPT + fit_results[10]*fit_results[10]*SMDR_GFermi_EXPT_UNC*SMDR_GFermi_EXPT_UNC/SMDR_GFermi_EXPT/SMDR_GFermi_EXPT );
		
		//double Delta_g = g_0 * ( TMath::Abs(fit_results[0])*SMDR_Mh_EXPT_UNC/SMDR_Mh_EXPT + TMath::Abs(fit_results[1])*SMDR_Mt_EXPT_UNC/SMDR_Mt_EXPT + TMath::Abs(fit_results[2])*SMDR_MZ_EXPT_UNC/SMDR_MZ_EXPT + TMath::Abs(fit_results[3])*SMDR_alphaS_MZ_EXPT_UNC/SMDR_alphaS_MZ_EXPT + TMath::Abs(fit_results[9])*SMDR_Delta_alpha_had_5_MZ_EXPT_UNC/SMDR_Delta_alpha_had_5_MZ_EXPT + TMath::Abs(fit_results[10])*SMDR_GFermi_EXPT_UNC/SMDR_GFermi_EXPT );
		
		cout << "\n++++++++++++++++++++++++++++++++++++++\nPropogation of errors: \nDelta g = " << Delta_g << "\nDelta x2 = " << g_0 * Delta_g / (8. * TMath::Pi()*TMath::Pi()) << endl; 
		
		outfile << "\n++++++++++++++++++++++++++++++++++++++\nPropogation of errors: \nDelta g = " << Delta_g << "\nDelta x2 = " << g_0 * Delta_g / (8. * TMath::Pi()*TMath::Pi()) << endl; 
		
		// Assessing
		
		cout << endl << "\n++++++++++++++++++++++++++++++++++++++\nASSESSING\n" << endl;
		tree_assess->SetBranchAddress("g",&g);
		
		double g_min = g_0;
		double g_max = g_0;
		double g_errmean = 0;
		double g_errmax = 0;
		
		double g_fromfit, g_err;
		
		int N = tree_assess->GetEntries();
		
		double onshell[7];
		
		for(int ientry = 0; ientry < N; ientry++){
			
			tree_assess->GetEntry(ientry);
			
			if(g < g_min) g_min = g;
			if(g > g_max) g_max = g;
			
			onshell[0] = Mh;
			onshell[1] = Mt;
			onshell[2] = MZ;
			onshell[3] = alphaS;
			onshell[4] = Delta;
			onshell[5] = mbmb;
			onshell[6] = GF;
			
			g_fromfit = fit_g(onshell, fit_results);
			
			g_err = TMath::Abs( g_fromfit - g );
			
			g_errmean += g_err / N;
			
			if (g_err > g_errmax) g_errmax = g_err;
			
			cout << "Processed event: " << ientry+1 << "/" << N << "\r";
			cout.flush();
		
		}
		
		cout << "\n\ng_mean = " << g_0 << "\ng_min = " << g_min << "\ng_max = " << g_max << "\nsigma = " << Delta_g << " \nParametric error: " << Delta_g / g_0 << "\n\ng_errmean = " << g_errmean << "\ng_errmax = " << g_errmax << "\nFit mean error: " << g_errmean / g_0 << "\nFit max error: " << g_errmax / g_0 << endl; 
		
		outfile << "\nParametric error: " << Delta_g / g_0 << "\ng_errmean = " << g_errmean << "\ng_errmax = " << g_errmax << "\nFit mean error: " << g_errmean / g_0 << "\nFit max error: " << g_errmax / g_0 << endl; 
		
		cout << "\n++++++++++++++++++++++++++++++++++++++\nFor table: \n{ " << fit_results[11];
		for (int k = 0; k < 11; k++) cout << ", " << fit_results[k];
		cout << ", " << Delta_g / g_0 << ", " << g_errmean / g_0 << ", " << g_errmax / g_0 << " }" << endl;
		
		outfile << "\nFor table: \n{ " << fit_results[11];
		for (int k = 0; k < 11; k++) outfile << ", " << fit_results[k];
		outfile << ", " << Delta_g / g_0 << ", " << g_errmean / g_0 << ", " << g_errmax / g_0 << " }\n++++++++++++++++++++++++++++++++++++++" << endl;
		
	}
	
	if(what==2)
	{
		
		tree->SetBranchAddress("gp",&gp);

		ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");			
		minimizer->SetMaxFunctionCalls(10000);
		minimizer->SetMaxIterations(10000);
		minimizer->SetTolerance(0.000001);
		minimizer->SetPrintLevel(1);

		ROOT::Math::Functor functor(&sum_squares_gp, 10); // <---!
		minimizer->SetFunction(functor);

		minimizer->SetVariable(0, "c h", 2.6e-7, 1e-8);
		minimizer->SetVariable(1, "c t", -2.609e-5, 1e-8);
		minimizer->SetVariable(2, "c Z", -4.7e-7, 1e-9);
		minimizer->SetVariable(3, "c S", 3.29e-6, 1e-8);
		minimizer->SetVariable(4, "c tt", -5e-8, 1e-10);
		minimizer->SetVariable(5, "c Zt", 1e-9, 1e-11);
		minimizer->SetVariable(6, "c ZZ", 1e-10, 1e-11);
		minimizer->SetVariable(7, "c a", 7.714e-5, 1e-8);
		minimizer->SetVariable(8, "c GF", 0.5, 0.001);			
		minimizer->SetVariable(9, "gp 0", 0.3586, 1e-5);
		
		minimizer->Minimize();
		const double* fit_results = minimizer->X();
		
		double gp_0 = fit_results[9];
		
		//Propogation of errors
		
		double Delta_gp = gp_0 * TMath::Sqrt( fit_results[0]*fit_results[0]*SMDR_Mh_EXPT_UNC*SMDR_Mh_EXPT_UNC/SMDR_Mh_EXPT/SMDR_Mh_EXPT + fit_results[1]*fit_results[1]*SMDR_Mt_EXPT_UNC*SMDR_Mt_EXPT_UNC/SMDR_Mt_EXPT/SMDR_Mt_EXPT + fit_results[2]*fit_results[2]*SMDR_MZ_EXPT_UNC*SMDR_MZ_EXPT_UNC/SMDR_MZ_EXPT/SMDR_MZ_EXPT + fit_results[3]*fit_results[3]*SMDR_alphaS_MZ_EXPT_UNC*SMDR_alphaS_MZ_EXPT_UNC/SMDR_alphaS_MZ_EXPT/SMDR_alphaS_MZ_EXPT + fit_results[7]*fit_results[7]*SMDR_Delta_alpha_had_5_MZ_EXPT_UNC*SMDR_Delta_alpha_had_5_MZ_EXPT_UNC/SMDR_Delta_alpha_had_5_MZ_EXPT/SMDR_Delta_alpha_had_5_MZ_EXPT + fit_results[8]*fit_results[8]*SMDR_GFermi_EXPT_UNC*SMDR_GFermi_EXPT_UNC/SMDR_GFermi_EXPT/SMDR_GFermi_EXPT );
		
		//double Delta_gp = gp_0 * ( TMath::Abs(fit_results[0])*SMDR_Mh_EXPT_UNC/SMDR_Mh_EXPT + TMath::Abs(fit_results[1])*SMDR_Mt_EXPT_UNC/SMDR_Mt_EXPT + TMath::Abs(fit_results[2])*SMDR_MZ_EXPT_UNC/SMDR_MZ_EXPT + TMath::Abs(fit_results[3])*SMDR_alphaS_MZ_EXPT_UNC/SMDR_alphaS_MZ_EXPT + TMath::Abs(fit_results[7])*SMDR_Delta_alpha_had_5_MZ_EXPT_UNC/SMDR_Delta_alpha_had_5_MZ_EXPT + TMath::Abs(fit_results[8])*SMDR_GFermi_EXPT_UNC/SMDR_GFermi_EXPT );
		
		cout << "\n++++++++++++++++++++++++++++++++++++++\nPropogation of errors: \nDelta gp = " << Delta_gp << "\nDelta x1 = " << 5. * gp_0 * Delta_gp / (24. * TMath::Pi()*TMath::Pi()) << endl;
		
		outfile << "\n++++++++++++++++++++++++++++++++++++++\nPropogation of errors: \nDelta gp = " << Delta_gp << "\nDelta x1 = " << 5. * gp_0 * Delta_gp / (24. * TMath::Pi()*TMath::Pi()) << endl; 
		
		// Assessing
		
		cout << endl << "\n++++++++++++++++++++++++++++++++++++++\nASSESSING\n" << endl;
		tree_assess->SetBranchAddress("gp",&gp);
		
		double gp_min = gp_0;
		double gp_max = gp_0;
		double gp_errmean = 0;
		double gp_errmax = 0;
		
		double gp_fromfit, gp_err;
		
		int N = tree_assess->GetEntries();
		
		double onshell[7];
		
		for(int ientry = 0; ientry < N; ientry++){
			
			tree_assess->GetEntry(ientry);
			
			if(gp < gp_min) gp_min = gp;
			if(gp > gp_max) gp_max = gp;
			
			onshell[0] = Mh;
			onshell[1] = Mt;
			onshell[2] = MZ;
			onshell[3] = alphaS;
			onshell[4] = Delta;
			onshell[5] = mbmb;
			onshell[6] = GF;
			
			gp_fromfit = fit_gp(onshell, fit_results);
			
			gp_err = TMath::Abs( gp_fromfit - gp );
			
			gp_errmean += gp_err / N;
			
			if (gp_err > gp_errmax) gp_errmax = gp_err;
			
			cout << "Processed event: " << ientry+1 << "/" << N << "\r";
			cout.flush();
		
		}
		
		cout << "\n\ngp_mean = " << gp_0 << "\ngp_min = " << gp_min << "\ngp_max = " << gp_max << "\nsigma = " << Delta_gp << " \nParametric error: " << Delta_gp / gp_0 << "\n\ngp_errmean = " << gp_errmean << "\ngp_errmax = " << gp_errmax << "\nFit mean error: " << gp_errmean / gp_0 << "\nFit max error: " << gp_errmax / gp_0 << endl; 
		
		outfile << "\nParametric error: " << Delta_gp / gp_0 << "\ngp_errmean = " << gp_errmean << "\ngp_errmax = " << gp_errmax << "\nFit mean error: " << gp_errmean / gp_0 << "\nFit max error: " << gp_errmax / gp_0 << endl; 
		
		cout << "\n++++++++++++++++++++++++++++++++++++++\nFor table: \n{ " << fit_results[9];
		for (int k = 0; k < 9; k++) cout << ", " << fit_results[k];
		cout << ", " << Delta_gp / gp_0 << ", " << gp_errmean / gp_0 << ", " << gp_errmax / gp_0 << " }" << endl;
		
		outfile << "\nFor table: \n{ " << fit_results[9];
		for (int k = 0; k < 9; k++) outfile << ", " << fit_results[k];
		outfile << ", " << Delta_gp / gp_0 << ", " << gp_errmean / gp_0 << ", " << gp_errmax / gp_0 << " }\n++++++++++++++++++++++++++++++++++++++" << endl;
		
	}
	
	if(what==3)
	{
		
		tree->SetBranchAddress("g3",&g3);

		ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");			
		minimizer->SetMaxFunctionCalls(10000);
		minimizer->SetMaxIterations(10000);
		minimizer->SetTolerance(0.000001);
		minimizer->SetPrintLevel(1);

		ROOT::Math::Functor functor(&sum_squares_g3, 9); // <---!
		minimizer->SetFunction(functor);

		minimizer->SetVariable(0, "c h", 2.8e-8, 1e-9);
		minimizer->SetVariable(1, "c t", -3.98e-5, 1e-7);
		minimizer->SetVariable(2, "c Z", 2.7e-9, 1e-10);
		minimizer->SetVariable(3, "c S", 3.7875e-3, 1e-7);
		minimizer->SetVariable(4, "c tt", 1.5e-8, 1e-9); //??
		minimizer->SetVariable(5, "c tS", 1e-7, 1e-8); //??
		minimizer->SetVariable(6, "c SS", -1.07e-5, 1e-7);
		minimizer->SetVariable(7, "c SSS", 5e-8, 1e-10);
		minimizer->SetVariable(8, "g3 0", 1.117, 1e-4);
		
		minimizer->Minimize();
		const double* fit_results = minimizer->X();
		
		double g3_0 = fit_results[8];
		
		//Propogation of errors
		
		double Delta_g3 = g3_0 * TMath::Sqrt( fit_results[0]*fit_results[0]*SMDR_Mh_EXPT_UNC*SMDR_Mh_EXPT_UNC/SMDR_Mh_EXPT/SMDR_Mh_EXPT + fit_results[1]*fit_results[1]*SMDR_Mt_EXPT_UNC*SMDR_Mt_EXPT_UNC/SMDR_Mt_EXPT/SMDR_Mt_EXPT + fit_results[2]*fit_results[2]*SMDR_MZ_EXPT_UNC*SMDR_MZ_EXPT_UNC/SMDR_MZ_EXPT/SMDR_MZ_EXPT + fit_results[3]*fit_results[3]*SMDR_alphaS_MZ_EXPT_UNC*SMDR_alphaS_MZ_EXPT_UNC/SMDR_alphaS_MZ_EXPT/SMDR_alphaS_MZ_EXPT );
		
		//double Delta_g3 = g3_0 * ( TMath::Abs(fit_results[0])*SMDR_Mh_EXPT_UNC/SMDR_Mh_EXPT + TMath::Abs(fit_results[1])*SMDR_Mt_EXPT_UNC/SMDR_Mt_EXPT + TMath::Abs(fit_results[2])*SMDR_MZ_EXPT_UNC/SMDR_MZ_EXPT + TMath::Abs(fit_results[3])*SMDR_alphaS_MZ_EXPT_UNC/SMDR_alphaS_MZ_EXPT );
		
		cout << "\n++++++++++++++++++++++++++++++++++++++\nPropogation of errors: \nDelta g3 = " << Delta_g3 << "\nDelta x3 = " << g3_0 * Delta_g3 / (8. * TMath::Pi()*TMath::Pi()) << endl; 
		
		outfile << "\n++++++++++++++++++++++++++++++++++++++\nPropogation of errors: \nDelta g3 = " << Delta_g3 << "\nDelta x3 = " << g3_0 * Delta_g3 / (8. * TMath::Pi()*TMath::Pi()) << endl; 
		
		// Assessing
		
		cout << endl << "\n++++++++++++++++++++++++++++++++++++++\nASSESSING\n" << endl;
		tree_assess->SetBranchAddress("g3",&g3);
		
		double g3_min = g3_0;
		double g3_max = g3_0;
		double g3_errmean = 0;
		double g3_errmax = 0;
		
		double g3_fromfit, g3_err;
		
		int N = tree_assess->GetEntries();
		
		double onshell[7];
		
		for(int ientry = 0; ientry < N; ientry++){
			
			tree_assess->GetEntry(ientry);
			
			if(g3 < g3_min) g3_min = g3;
			if(g3 > g3_max) g3_max = g3;
			
			onshell[0] = Mh;
			onshell[1] = Mt;
			onshell[2] = MZ;
			onshell[3] = alphaS;
			onshell[4] = Delta;
			onshell[5] = mbmb;
			onshell[6] = GF;
			
			g3_fromfit = fit_g3(onshell, fit_results);
			
			g3_err = TMath::Abs( g3_fromfit - g3 );
			
			g3_errmean += g3_err / N;
			
			if (g3_err > g3_errmax) g3_errmax = g3_err;
			
			cout << "Processed event: " << ientry+1 << "/" << N << "\r";
			cout.flush();
		
		}
		
		cout << "\n\ng3_mean = " << g3_0 << "\ng3_min = " << g3_min << "\ng3_max = " << g3_max << "\nsigma = " << Delta_g3 << " \nParametric error: " << Delta_g3 / g3_0 << "\n\ng3_errmean = " << g3_errmean << "\ng3_errmax = " << g3_errmax << "\nFit mean error: " << g3_errmean / g3_0 << "\nFit max error: " << g3_errmax / g3_0 << endl; 
		
		outfile << "\nParametric error: " << Delta_g3 / g3_0 << "\ng3_errmean = " << g3_errmean << "\ng3_errmax = " << g3_errmax << "\nFit mean error: " << g3_errmean / g3_0 << "\nFit max error: " << g3_errmax / g3_0 << endl; 
		
		cout << "\n++++++++++++++++++++++++++++++++++++++\nFor table: \n{ " << fit_results[8];
		for (int k = 0; k < 8; k++) cout << ", " << fit_results[k];
		cout << ", " << Delta_g3 / g3_0 << ", " << g3_errmean / g3_0 << ", " << g3_errmax / g3_0 << " }" << endl;
		
		outfile << "\nFor table: \n{ " << fit_results[8];
		for (int k = 0; k < 8; k++) outfile << ", " << fit_results[k];
		outfile << ", " << Delta_g3 / g3_0 << ", " << g3_errmean / g3_0 << ", " << g3_errmax / g3_0 << " }\n++++++++++++++++++++++++++++++++++++++" << endl;
		
	}
	
	if(what==4)
	{
		
		tree->SetBranchAddress("yt",&yt);

		ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");			
		minimizer->SetMaxFunctionCalls(10000);
		minimizer->SetMaxIterations(10000);
		minimizer->SetTolerance(0.000001);
		minimizer->SetPrintLevel(1);

		ROOT::Math::Functor functor(&sum_squares_yt, 10); // <---!
		minimizer->SetFunction(functor);

		minimizer->SetVariable(0, "c h", -2.36e-6, 1e-8);
		minimizer->SetVariable(1, "c t", 6.352e-3, 1e-6);
		minimizer->SetVariable(2, "c Z", -3.5e-7, 1e-8);
		minimizer->SetVariable(3, "c S", -7.76e-4, 1e-6);
		minimizer->SetVariable(4, "c tt", 1.5e-9, 1e-10);
		minimizer->SetVariable(5, "c tS", 1e-8, 1e-9); 
		minimizer->SetVariable(6, "c SS", 3e-7, 1e-8);
		minimizer->SetVariable(7, "c ht", -5e-8, 1e-10);
		minimizer->SetVariable(8, "c a", 2.2e-8, 1e-9);	
		minimizer->SetVariable(9, "yt 0", 0.9496, 1e-5);
		
		minimizer->Minimize();
		const double* fit_results = minimizer->X();
		
		double yt_0 = fit_results[9];
		
		//Propogation of errors
		
		double Delta_yt = yt_0 * TMath::Sqrt( fit_results[0]*fit_results[0]*SMDR_Mh_EXPT_UNC*SMDR_Mh_EXPT_UNC/SMDR_Mh_EXPT/SMDR_Mh_EXPT + fit_results[1]*fit_results[1]*SMDR_Mt_EXPT_UNC*SMDR_Mt_EXPT_UNC/SMDR_Mt_EXPT/SMDR_Mt_EXPT + fit_results[2]*fit_results[2]*SMDR_MZ_EXPT_UNC*SMDR_MZ_EXPT_UNC/SMDR_MZ_EXPT/SMDR_MZ_EXPT + fit_results[3]*fit_results[3]*SMDR_alphaS_MZ_EXPT_UNC*SMDR_alphaS_MZ_EXPT_UNC/SMDR_alphaS_MZ_EXPT/SMDR_alphaS_MZ_EXPT + fit_results[8]*fit_results[8]*SMDR_Delta_alpha_had_5_MZ_EXPT_UNC*SMDR_Delta_alpha_had_5_MZ_EXPT_UNC/SMDR_Delta_alpha_had_5_MZ_EXPT/SMDR_Delta_alpha_had_5_MZ_EXPT );
		
		//double Delta_yt = yt_0 * ( TMath::Abs(fit_results[0])*SMDR_Mh_EXPT_UNC/SMDR_Mh_EXPT + TMath::Abs(fit_results[1])*SMDR_Mt_EXPT_UNC/SMDR_Mt_EXPT + TMath::Abs(fit_results[2])*SMDR_MZ_EXPT_UNC/SMDR_MZ_EXPT + TMath::Abs(fit_results[3])*SMDR_alphaS_MZ_EXPT_UNC/SMDR_alphaS_MZ_EXPT + TMath::Abs(fit_results[8])*SMDR_Delta_alpha_had_5_MZ_EXPT_UNC/SMDR_Delta_alpha_had_5_MZ_EXPT );
		
		cout << "\n++++++++++++++++++++++++++++++++++++++\nPropogation of errors: \nDelta yt = " << Delta_yt << "\nDelta y1 = " << yt_0 * Delta_yt / (8. * TMath::Pi()*TMath::Pi()) << endl; 
		
		outfile << "\n++++++++++++++++++++++++++++++++++++++\nPropogation of errors: \nDelta yt = " << Delta_yt << "\nDelta y1 = " << yt_0 * Delta_yt / (8. * TMath::Pi()*TMath::Pi()) << endl; 
		
		// Assessing
		
		cout << endl << "\n++++++++++++++++++++++++++++++++++++++\nASSESSING\n" << endl;
		tree_assess->SetBranchAddress("yt",&yt);
		
		double yt_min = yt_0;
		double yt_max = yt_0;
		double yt_errmean = 0;
		double yt_errmax = 0;
		
		double yt_fromfit, yt_err;
		
		int N = tree_assess->GetEntries();
		
		double onshell[7];
		
		for(int ientry = 0; ientry < N; ientry++){
			
			tree_assess->GetEntry(ientry);
			
			if(yt < yt_min) yt_min = yt;
			if(yt > yt_max) yt_max = yt;
			
			onshell[0] = Mh;
			onshell[1] = Mt;
			onshell[2] = MZ;
			onshell[3] = alphaS;
			onshell[4] = Delta;
			onshell[5] = mbmb;
			onshell[6] = GF;
			
			yt_fromfit = fit_yt(onshell, fit_results);
			
			yt_err = TMath::Abs( yt_fromfit - yt );
			
			yt_errmean += yt_err / N;
			
			if (yt_err > yt_errmax) yt_errmax = yt_err;
			
			cout << "Processed event: " << ientry+1 << "/" << N << "\r";
			cout.flush();
		
		}
		
		cout << "\n\nyt_mean = " << yt_0 << "\nyt_min = " << yt_min << "\nyt_max = " << yt_max << "\nsigma = " << Delta_yt << " \nParametric error: " << Delta_yt / yt_0 << "\n\nyt_errmean = " << yt_errmean << "\nyt_errmax = " << yt_errmax << "\nFit mean error: " << yt_errmean / yt_0 << "\nFit max error: " << yt_errmax / yt_0 << endl; 
		
		outfile << "\nParametric error: " << Delta_yt / yt_0 << "\nyt_errmean = " << yt_errmean << "\nyt_errmax = " << yt_errmax << "\nFit mean error: " << yt_errmean / yt_0 << "\nFit max error: " << yt_errmax / yt_0 << endl; 
		
		cout << "\n++++++++++++++++++++++++++++++++++++++\nFor table: \n{ " << fit_results[9];
		for (int k = 0; k < 9; k++) cout << ", " << fit_results[k];
		cout << ", " << Delta_yt / yt_0 << ", " << yt_errmean / yt_0 << ", " << yt_errmax / yt_0 << " }" << endl;
		
		outfile << "\nFor table: \n{ " << fit_results[9];
		for (int k = 0; k < 9; k++) outfile << ", " << fit_results[k];
		outfile << ", " << Delta_yt / yt_0 << ", " << yt_errmean / yt_0 << ", " << yt_errmax / yt_0 << " }\n++++++++++++++++++++++++++++++++++++++" << endl;
		
	}
	
	if(what==5)
	{
		
		tree->SetBranchAddress("yb",&yb);

		ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");			
		minimizer->SetMaxFunctionCalls(10000);
		minimizer->SetMaxIterations(10000);
		minimizer->SetTolerance(0.000001);
		minimizer->SetPrintLevel(1);

		ROOT::Math::Functor functor(&sum_squares_yb, 11); // <---!
		minimizer->SetFunction(functor);
		
		minimizer->SetVariable(0, "c b", 1.185, 1e-3);
		minimizer->SetVariable(1, "c h", -2e-8, 1e-9);
		minimizer->SetVariable(2, "c t", -2.4e-5, 1e-6);
		minimizer->SetVariable(3, "c Z", -2e-8, 1e-9);
		minimizer->SetVariable(4, "c S", -6.125e-3, 1e-6);
		minimizer->SetVariable(5, "c bb", 0.075, 1e-3);
		minimizer->SetVariable(6, "c bS", -3.3e-3, 1e-4);
		minimizer->SetVariable(7, "c SS", -2.1e-5, 1e-6);
		minimizer->SetVariable(8, "c tS", -2e-8, 1e-9);
		minimizer->SetVariable(9, "c SSS", -2e-8, 1e-9);
		minimizer->SetVariable(10, "yb 0", 0.01684, 1e-6);
		
		minimizer->Minimize();
		const double* fit_results = minimizer->X();
		
		double yb_0 = fit_results[10];
		
		//Propogation of errors
		
		double Delta_yb = yb_0 * TMath::Sqrt( fit_results[0]*fit_results[0]*SMDR_mbmb_EXPT_UNC_hi*SMDR_mbmb_EXPT_UNC_hi/SMDR_mbmb_EXPT/SMDR_mbmb_EXPT + fit_results[1]*fit_results[1]*SMDR_Mh_EXPT_UNC*SMDR_Mh_EXPT_UNC/SMDR_Mh_EXPT/SMDR_Mh_EXPT + fit_results[2]*fit_results[2]*SMDR_Mt_EXPT_UNC*SMDR_Mt_EXPT_UNC/SMDR_Mt_EXPT/SMDR_Mt_EXPT + fit_results[3]*fit_results[3]*SMDR_MZ_EXPT_UNC*SMDR_MZ_EXPT_UNC/SMDR_MZ_EXPT/SMDR_MZ_EXPT + fit_results[4]*fit_results[4]*SMDR_alphaS_MZ_EXPT_UNC*SMDR_alphaS_MZ_EXPT_UNC/SMDR_alphaS_MZ_EXPT/SMDR_alphaS_MZ_EXPT );
		
		//double Delta_yb = yb_0 * ( TMath::Abs(fit_results[0])*SMDR_mbmb_EXPT_UNC_hi/SMDR_mbmb_EXPT + TMath::Abs(fit_results[1])*SMDR_Mh_EXPT_UNC/SMDR_Mh_EXPT + TMath::Abs(fit_results[2])*SMDR_Mt_EXPT_UNC/SMDR_Mt_EXPT + TMath::Abs(fit_results[3])*SMDR_MZ_EXPT_UNC/SMDR_MZ_EXPT + TMath::Abs(fit_results[4])*SMDR_alphaS_MZ_EXPT_UNC/SMDR_alphaS_MZ_EXPT );
		
		cout << "\n++++++++++++++++++++++++++++++++++++++\nPropogation of errors: \nDelta yb = " << Delta_yb << "\nDelta y2 = " << yb_0 * Delta_yb / (8. * TMath::Pi()*TMath::Pi()) << endl; 
		
		outfile << "\n++++++++++++++++++++++++++++++++++++++\nPropogation of errors: \nDelta yb = " << Delta_yb << "\nDelta y2 = " << yb_0 * Delta_yb / (8. * TMath::Pi()*TMath::Pi()) << endl; 
		
		// Assessing
		
		cout << endl << "\n++++++++++++++++++++++++++++++++++++++\nASSESSING\n" << endl;
		tree_assess->SetBranchAddress("yb",&yb);
		
		double yb_min = yb_0;
		double yb_max = yb_0;
		double yb_errmean = 0;
		double yb_errmax = 0;
		
		double yb_fromfit, yb_err;
		
		int N = tree_assess->GetEntries();
		
		double onshell[7];
		
		for(int ientry = 0; ientry < N; ientry++){
			
			tree_assess->GetEntry(ientry);
			
			if(yb < yb_min) yb_min = yb;
			if(yb > yb_max) yb_max = yb;
			
			onshell[0] = Mh;
			onshell[1] = Mt;
			onshell[2] = MZ;
			onshell[3] = alphaS;
			onshell[4] = Delta;
			onshell[5] = mbmb;
			onshell[6] = GF;
			
			yb_fromfit = fit_yb(onshell, fit_results);
			
			yb_err = TMath::Abs( yb_fromfit - yb );
			
			yb_errmean += yb_err / N;
			
			if (yb_err > yb_errmax) yb_errmax = yb_err;
			
			cout << "Processed event: " << ientry+1 << "/" << N << "\r";
			cout.flush();
		
		}
		
		cout << "\n\nyb_mean = " << yb_0 << "\nyb_min = " << yb_min << "\nyb_max = " << yb_max << "\nsigma = " << Delta_yb << " \nParametric error: " << Delta_yb / yb_0 << "\n\nyb_errmean = " << yb_errmean << "\nyb_errmax = " << yb_errmax << "\nFit mean error: " << yb_errmean / yb_0 << "\nFit max error: " << yb_errmax / yb_0 << endl; 
		
		outfile << "\nParametric error: " << Delta_yb / yb_0 << "\nyb_errmean = " << yb_errmean << "\nyb_errmax = " << yb_errmax << "\nFit mean error: " << yb_errmean / yb_0 << "\nFit max error: " << yb_errmax / yb_0 << endl; 
		
		cout << "\n++++++++++++++++++++++++++++++++++++++\nFor table: \n{ " << fit_results[10];
		for (int k = 0; k < 10; k++) cout << ", " << fit_results[k];
		cout << ", " << Delta_yb / yb_0 << ", " << yb_errmean / yb_0 << ", " << yb_errmax / yb_0 << " }" << endl;
		
		outfile << "\nFor table: \n{ " << fit_results[10];
		for (int k = 0; k < 10; k++) outfile << ", " << fit_results[k];
		outfile << ", " << Delta_yb / yb_0 << ", " << yb_errmean / yb_0 << ", " << yb_errmax / yb_0 << " }\n++++++++++++++++++++++++++++++++++++++" << endl;
		
	}
	
	outfile.close();	
	
	return 0;

}
