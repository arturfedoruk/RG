// basically generate_net but with loop_Fit_Inputs_original()
// to launch the program: 
// g++ `root-config --cflags` generate_random_orig.cpp `root-config --libs` -lm -lsmdr -ltsil -l3vil

#include "smdr.h"
#include "iostream"
#include "fstream"
#include "cstdio"
#include "string"
#include "TROOT.h"
#include "TH1D.h"
#include "TError.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
using namespace std;

#include "../loop_Fit_Inputs_original.cpp"
#include "../loop_configs.cpp"

int main(){
	
	/*
	int seed;
	ifstream lastSeed("lastSeed.txt");
	if ( lastSeed.is_open() ) {
		lastSeed >> seed; 
	}
	lastSeed.close();
	remove("lastSeed.txt");
	ofstream lastSeed_new("lastSeed.txt");
	lastSeed_new << seed+1 << endl;
	lastSeed_new.close();
	*/
	TRandom3 r(0);
	int seed = int( r.Rndm() * 2147483640 );
	
	
	TFile *file = new TFile("random_data_orig.root", "UPDATE");
	TTree *tree;
	
	int input_config, N;
	
	cout << "\n\nEnter number of events to be generated: ";
	cin >> N;
	
	bool alright = true;
	
	if (N < 1) alright = false;
	if (!alright) cout << "ERROR: number of events must be >0 !!!" << endl;
	
	float config[9];
	if(alright){
	
	cout << "Enter loop configuration: ";
	cin >> input_config;
	
	if (input_config == 111111){
		for (int k = 0; k < 9; k++) config[k] = config_111111[k];
		tree = (TTree*)gDirectory->Get("tree_111111");
	} else if (input_config == 222222){
		for (int k = 0; k < 9; k++) config[k] = config_222222[k];
		tree = (TTree*)gDirectory->Get("tree_222222");
	} else if (input_config == 333333){
		for (int k = 0; k < 9; k++) config[k] = config_333333[k];
		tree = (TTree*)gDirectory->Get("tree_333333");
	} else if (input_config == 333221){
		for (int k = 0; k < 9; k++) config[k] = config_333221[k];
		tree = (TTree*)gDirectory->Get("tree_333221");
	} else if (input_config == 444332){
		for (int k = 0; k < 9; k++) config[k] = config_444332[k];
		tree = (TTree*)gDirectory->Get("tree_444332");
	} else if (input_config == 444333){
		for (int k = 0; k < 9; k++) config[k] = config_444333[k];
		tree = (TTree*)gDirectory->Get("tree_444333");
	} else {
		alright = false;
		cout << "ERROR: such a configuration is not supported, supported ones: 111111, 222222, 333333, 333221, 444332, 444333" << endl;
	}
	
	}
	
	if (alright){
	cout << "\nGeneration will be performed using TRandom3 with seed " << seed << endl;
	
	double alphaS_MZ, alpha, GFermi, MZ, Mh, Mt, mbmb, Delta_alpha, g, gp, g3, yt, yb, lambda, x1, x2, x3, y1, y2, z;
	tree->SetBranchAddress("seed",&seed);
	tree->SetBranchAddress("alphaS_MZ",&alphaS_MZ);
	tree->SetBranchAddress("alpha",&alpha);
	tree->SetBranchAddress("GFermi",&GFermi);
	tree->SetBranchAddress("MZ",&MZ);
	tree->SetBranchAddress("Mh",&Mh);
	tree->SetBranchAddress("Mt",&Mt);
	tree->SetBranchAddress("mbmb",&mbmb);
	tree->SetBranchAddress("Delta_alpha",&Delta_alpha);
	tree->SetBranchAddress("g",&g);
	tree->SetBranchAddress("gp",&gp);
	tree->SetBranchAddress("g3",&g3);
	tree->SetBranchAddress("yt",&yt);
	tree->SetBranchAddress("yb",&yb);
	tree->SetBranchAddress("lambda",&lambda);
	tree->SetBranchAddress("x1",&x1);
	tree->SetBranchAddress("x2",&x2);
	tree->SetBranchAddress("x3",&x3);
	tree->SetBranchAddress("y1",&y1);
	tree->SetBranchAddress("y2",&y2);
	tree->SetBranchAddress("z",&z);
	
	TRandom3 random_generator(seed);
	double ERROR_TOLERANCE = 1e-6;
	double pi = 3.14159265359;
	
	for (int ientry = 0; ientry < N; ientry++){
		SMDR_Q_in = 200;
	
		alphaS_MZ = random_generator.Gaus(SMDR_alphaS_MZ_EXPT, SMDR_alphaS_MZ_EXPT_UNC);
		alpha = random_generator.Gaus(SMDR_alpha_EXPT, SMDR_alpha_EXPT_UNC);
		GFermi = random_generator.Gaus(SMDR_GFermi_EXPT, SMDR_GFermi_EXPT_UNC);
		MZ = random_generator.Gaus(SMDR_MZ_EXPT, SMDR_MZ_EXPT_UNC);
		Mh = random_generator.Gaus(SMDR_Mh_EXPT, SMDR_Mh_EXPT_UNC);
		Mt = random_generator.Gaus(SMDR_Mt_EXPT, SMDR_Mt_EXPT_UNC);
		mbmb = random_generator.Gaus(SMDR_mbmb_EXPT, SMDR_mbmb_EXPT_UNC_hi);
		Delta_alpha = random_generator.Gaus(SMDR_Delta_alpha_had_5_MZ_EXPT, SMDR_Delta_alpha_had_5_MZ_EXPT_UNC);
		
		loop_Fit_Inputs_original (SMDR_Q_in, 
				alphaS_MZ, 
				alpha, GFermi, 
				MZ, Mh, Mt, 
				mbmb, Delta_alpha,
		           	ERROR_TOLERANCE,
		           	config);
		           	
	      	g = SMDR_g_in;
		gp = SMDR_gp_in;
		g3 = SMDR_g3_in;
		yt = SMDR_yt_in;
		yb = SMDR_yb_in;
		lambda = SMDR_lambda_in;

		x1 = SMDR_gp_in*SMDR_gp_in*5 / (48*pi*pi);
		x2 = SMDR_g_in*SMDR_g_in / (16*pi*pi);
		x3 = SMDR_g3_in*SMDR_g3_in / (16*pi*pi);
		y1 = SMDR_yt_in*SMDR_yt_in / (16*pi*pi);
		y2 = SMDR_yb_in*SMDR_yb_in / (16*pi*pi);
		z = SMDR_lambda_in / (16*pi*pi);
		
		tree->Fill();
		cout << "Processed event: " << ientry+1 << "/" << N << "\r";
		cout.flush();
	
	}
	
	tree->Write(0,TObject::kWriteDelete,0);
	cout << "\nSuccessfully processed " << N << " events" << endl;
	
	}

	file->Close("R");
	delete file;

	return 0;

}
