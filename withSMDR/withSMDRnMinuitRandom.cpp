// to launch the program: 
// g++ `root-config --cflags` withSMDRnMinuitRandom.cpp `root-config --libs` -lm -lsmdr -ltsil -l3vil

#include "smdr.h"
#include "iostream"
#include "fstream"
#include "string"
#include "TROOT.h"
#include "TH1D.h"
#include "TError.h"
#include "TRandom3.h"
#include "TCanvas.h"
using namespace std;
#define ZEROSAFE(a) (((a) > (SMDR_TOL)) ? (a) : (SMDR_TOL)) //idk wht's that

#include "my_Fit_Inputs.cpp"

int main(){
	
	double pi = 3.14159265359;
	
	int SEED = 1;
	
	TRandom3 random_generator(SEED);
	
	TH1D *h_x1 = new TH1D("h_x1", "random distribution of x1", 100, 0.0012, 0.0015);
	TH1D *h_x2 = new TH1D("h_x2", "random distribution of x2", 100, 0.0025, 0.0028);
	TH1D *h_x3 = new TH1D("h_x3", "random distribution of x3", 100, 0.0084, 0.0088);
	TH1D *h_y1 = new TH1D("h_y1", "random distribution of y1", 100, 0.0060, 0.0064);
	TH1D *h_y2 = new TH1D("h_y2", "random distribution of y2", 100, 1.6e-6, 1.9e-6);
	TH1D *h_z = new TH1D("h_z", "random distribution of z", 100, 0.00079, 0.00082);
	TH1D *h_check = new TH1D("h_Mt", "random distribution of Mt", 100, SMDR_Mt_EXPT - 0.7, SMDR_Mt_EXPT + 0.7);
	
	SMDR_Q_in = 173.22; // mass of top-quark, we work here
	float config_111111[9] = {0, 0, 0, 0, 0, 0, 1, 1} ;
	float config_222222[9] = {1, 1, 1, 1, 1, 1, 1, 1} ;
	SMDR_REAL Mt_tmp;
	
	int N = 1000;
	double ERROR_TOLERANCE = 1.e-6; // idk
	
	cout << "number of ivents: " << N << endl;
	
	for (int i = 0; i < N; i++){
		SMDR_Q_in = 173.22;
		Mt_tmp = random_generator.Gaus(SMDR_Mt_EXPT, SMDR_Mt_EXPT_UNC);
		
		my_Fit_Inputs (SMDR_Q_in,
		           random_generator.Gaus(SMDR_alphaS_MZ_EXPT, SMDR_alphaS_MZ_EXPT_UNC),
		           random_generator.Gaus(SMDR_alpha_EXPT, SMDR_alpha_EXPT_UNC),
		           random_generator.Gaus(SMDR_GFermi_EXPT, SMDR_GFermi_EXPT_UNC),
		           random_generator.Gaus(SMDR_MZ_EXPT, SMDR_MZ_EXPT_UNC),
		           random_generator.Gaus(SMDR_Mh_EXPT, SMDR_Mh_EXPT_UNC),
		           Mt_tmp,
		           random_generator.Gaus(SMDR_mbmb_EXPT, SMDR_mbmb_EXPT_UNC_hi),
		           random_generator.Gaus(SMDR_Delta_alpha_had_5_MZ_EXPT, SMDR_Delta_alpha_had_5_MZ_EXPT),
		           ERROR_TOLERANCE,
		           config_111111);
		           
	         h_check->Fill(Mt_tmp);
	         h_x1->Fill(SMDR_gp_in*SMDR_gp_in*5 / (48*pi*pi));
	         h_x2->Fill(SMDR_g_in*SMDR_g_in / (16*pi*pi));
	         h_x3->Fill(SMDR_g3_in*SMDR_g3_in / (16*pi*pi));
	         h_y1->Fill(SMDR_yt_in*SMDR_yt_in / (16*pi*pi));
	         h_y2->Fill(SMDR_yb_in*SMDR_yb_in / (16*pi*pi));
	         h_z->Fill(SMDR_lambda_in / (16*pi*pi));
	         
	         cout << "Processed " << i << "/" << N << endl;
	}
	
	TCanvas *c = new TCanvas(" "," ");
	h_check->Draw();
	h_check->Fit("gaus");
	c->SaveAs("hists/h_check.pdf");
	h_x1->Draw();
	h_x1->Fit("gaus");
	c->SaveAs("hists/h_x1.pdf");
	h_x2->Draw();
	h_x2->Fit("gaus");
	c->SaveAs("hists/h_x2.pdf");
	h_x3->Draw();
	h_x3->Fit("gaus");
	c->SaveAs("hists/h_x3.pdf");
	h_y1->Draw();
	h_y1->Fit("gaus");
	c->SaveAs("hists/h_y1.pdf");
	h_y2->Draw();
	h_y2->Fit("gaus");
	c->SaveAs("hists/h_y2.pdf");
	h_z->Draw();
	h_z->Fit("gaus");
	c->SaveAs("hists/h_z.pdf");
         
         cout << "sigmas: \n x1: " << h_x1->GetStdDev() << "\n x2: " << h_x2->GetStdDev() << "\n x3: " << h_x3->GetStdDev() << "\n y1: " << h_y1->GetStdDev() << "\n y2: " << h_y2->GetStdDev() << "\n z: " << h_z->GetStdDev() << endl;
	
	return 0;
}


