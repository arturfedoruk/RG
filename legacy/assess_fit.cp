// old program for assessing fits, now incorporated in fit.cpp
// to launch the program: 
// g++ `root-config --cflags` withSMDRassessfit.cpp `root-config --libs` -lm -lsmdr -ltsil -l3vil

#include "smdr.h"
#include "iostream"
#include "fstream"
#include "cstdio"
#include "string"
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TF1.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom3.h"
#include "TError.h"
using namespace std;

TFile *file = new TFile("net_data.root");
TTree *tree = (TTree*)file->Get("tree_222222_9,_5_5_5_5_5");

double fit_lambda(double Mh, double Mt, double MZ, double alphaS){
	double deltah = (Mh - SMDR_Mh_EXPT)/0.1;
	double deltat = (Mt - SMDR_Mt_EXPT);
	double deltaZ = (MZ - SMDR_MZ_EXPT)/0.001;
	double deltaS = (alphaS - SMDR_alphaS_MZ_EXPT)*1000;
	
	//return 0.128145 + 0.000205183 * deltah + 0.000119129 * deltat - 4.60135e-08 * deltaZ - 6.49278e-06 * deltaS; 
	return 0.127002 + 0.000205183 * deltah + 0.000105715 * deltat - 4.60135e-08 * deltaZ - 5.67882e-06 * deltaS + 3.45336e-06 * deltat*deltat - 3.38081e-07 * deltat*deltaS + 8.36871e-08 * deltah*deltah - 9.72156e-08 * deltah*deltat + 1.75564e-08 * deltaS*deltaS + 5.02476e-09 * deltah*deltaS + 3.04188e-08 * deltat*deltat*deltat - 3.29597e-09 * deltat*deltat*deltaS;
}

double fit_lambda2(double Mh, double Mt, double MZ, double alphaS, double GF, double Delta, double mbmb){
	double deltah = (Mh - SMDR_Mh_EXPT)/0.1;
	double deltat = (Mt - SMDR_Mt_EXPT);
	double deltaZ = (MZ - SMDR_MZ_EXPT)/0.001;
	double deltaS = (alphaS - SMDR_alphaS_MZ_EXPT)*1000;
	double deltaa = (Delta - SMDR_Delta_alpha_had_5_MZ_EXPT)*1e4;
	double Deltab = mbmb / SMDR_mbmb_EXPT - 1;
	double DeltaGF = GF / SMDR_GFermi_EXPT - 1;
	
	return 0.127001 * (1 + 0.00161514 * deltah + 0.000829787 * deltat - 3.61216e-07 * deltaZ - 4.43841e-05 * deltaS + 2.72037e-05 * deltat*deltat - 2.64946e-06 * deltat*deltaS + 6.58524e-07 * deltah*deltah - 7.66228e-07 * deltah*deltat + 1.38814e-07 * deltaS*deltaS + 3.93743e-08 * deltah*deltaS + 2.39648e-07 * deltat*deltat*deltat - 2.58566e-08 * deltat*deltat*deltaS + 3.28211e-07 * deltaa + 4.03967e-05 * Deltab + 0.997437 * DeltaGF);
}

double fit_lambda3(double Mh, double Mt, double MZ, double alphaS, double Delta){
	double deltah = (Mh - SMDR_Mh_EXPT)/0.1;
	double deltat = (Mt - SMDR_Mt_EXPT);
	double deltaZ = (MZ - SMDR_MZ_EXPT)/0.001;
	double deltaS = (alphaS - SMDR_alphaS_MZ_EXPT)*1000;
	double deltaa = (Delta - SMDR_Delta_alpha_had_5_MZ_EXPT)*1e4;
	
	return 0.127002 + 0.00161559 * deltah + 0.000938012 * deltat - 3.62305e-07 * deltaZ - 5.11236e-06 * deltaS + 2.71914e-05 * deltat*deltat - 2.66202e-06 * deltat*deltaS + 6.58943e-07 * deltah*deltah - 7.65465e-07 * deltah*deltat + 1.38237e-07 * deltaS*deltaS + 3.95644e-08 * deltah*deltaS + 3.4e-07 * deltaa;
}

int main(){

	double lambda_min = 0.127;
	double lambda_max = 0.127;
	double lambda_mean = 0;
	double lambda_errmean = 0;
	double lambda_errmax = 0;
	
	double alphaS, Mt, Mh, MZ, Delta, mbmb, GF, lambda, lambda_fromfit, lambda_err;
	
	tree->SetBranchAddress("alphaS_MZ",&alphaS);
	tree->SetBranchAddress("Mt",&Mt);
	tree->SetBranchAddress("Mh",&Mh);
	tree->SetBranchAddress("MZ",&MZ);
	
	tree->SetBranchAddress("Delta_alpha",&Delta);
	tree->SetBranchAddress("mbmb",&mbmb);
	tree->SetBranchAddress("GFermi",&GF);
	
	tree->SetBranchAddress("lambda",&lambda);
	
	int N = tree->GetEntries();
	
	for(int ientry = 0; ientry < N; ientry++){
		
		tree->GetEntry(ientry);
		
		lambda_mean += lambda / N;
		
		if(lambda < lambda_min) lambda_min = lambda;
		if(lambda > lambda_max) lambda_max = lambda;
		
		//lambda_fromfit = fit_lambda(Mh, Mt, MZ, alphaS);
		lambda_fromfit = fit_lambda2(Mh, Mt, MZ, alphaS, GF, Delta, mbmb);
		//lambda_fromfit = fit_lambda3(Mh, Mt, MZ, alphaS, Delta);
		lambda_err = TMath::Abs( lambda_fromfit - lambda );
		
		lambda_errmean += lambda_err / N;
		
		if (lambda_err > lambda_errmax) lambda_errmax = lambda_err;
		
		cout << "Processed event: " << ientry+1 << "/" << N << "\r";
		cout.flush();
	
	}
	
	cout << "\n\nlambda_mean = " << lambda_mean << "\nlambda_min = " << lambda_min << "\nlambda_max = " << lambda_max << "\nsigma = " << 0.000356 << " \nExperimental error: " << 0.000356 / lambda_mean << "\n\nlambda_errmean = " << lambda_errmean << "\nlambda_errmax = " << lambda_errmax << "\nFit max error: " << lambda_errmax / lambda_mean << endl;

	return 0;

}
