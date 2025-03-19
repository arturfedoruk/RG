// to launch the program: 
// g++ `root-config --cflags` withSMDRnROOTgeneratenetOrig.cpp `root-config --libs` -lm -lsmdr -ltsil -l3vil

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
#define ZEROSAFE(a) (((a) > (SMDR_TOL)) ? (a) : (SMDR_TOL)) //idk wht's that

#include "my_Fit_Inputs_original.cpp"

int main(){

	TFile *file = new TFile("net_data_orig.root","update");
	
	int input_config;	
	cout << "\nEnter loop configuration: ";
	cin >> input_config;
		
	float config_111111[9] = {0, 0, 0, 0, 0, 0, 1, 1} ; // for QCDQED_at_MZ & mbmb loop 0 doesn't exist
	float config_222222[9] = {1, 1, 1, 1, 1, 1, 1, 1} ;
	float config_333333[9] = {2, 2, 2, 2, 2, 2, 2, 2} ;
	float config_333221[9] = {2, 2, 0, 2, 2, 0, 2, 1} ;
	float config_444332[9] = {3, 2, 1, 2.5, 2.5, 1, 3, 2} ; 
	float config_444333[9] = {3, 2, 2, 2.5, 2.5, 2, 3, 2} ; // hfor MZ & MW loop 3 doesn't exist
	float config[9];
	
	TTree* tree;
	
	bool alright = true;
		
	if (input_config == 111111){
		for (int k = 0; k < 9; k++) config[k] = config_111111[k];
		tree = new TTree("tree_111111","tree_111111");
	} else if (input_config == 222222){
		for (int k = 0; k < 9; k++) config[k] = config_222222[k];
		tree = new TTree("tree_222222","tree_222222");
	} else if (input_config == 333333){
		for (int k = 0; k < 9; k++) config[k] = config_333333[k];
		tree = new TTree("tree_333333","tree_333333");
	} else if (input_config == 333221){
		for (int k = 0; k < 9; k++) config[k] = config_333221[k];
		tree = new TTree("tree_333221","tree_333221");
	} else if (input_config == 444332){
		for (int k = 0; k < 9; k++) config[k] = config_444332[k];
		tree = new TTree("tree_444332","tree_444332");
	} else if (input_config == 444333){
		for (int k = 0; k < 9; k++) config[k] = config_444333[k];
		tree = new TTree("tree_444333","tree_444333");
	} else {
		alright = false;
		cout << "ERROR: such a configuration is not supported, supported ones: 111111, 222222, 333333, 333221, 444332, 444333" << endl;
	}
	
	if(alright){
		
	double alphaS_MZ, alpha, GFermi, MZ, Mh, Mt, mbmb, Delta_alpha, g, gp, g3, yt, yb, lambda, x1, x2, x3, y1, y2, z;
	
	double ERROR_TOLERANCE = 1e-6;
	double pi = 3.14159265359;
	
	SMDR_Q_in = 200;
	
	cout << "Creating tree_" << input_config << endl;
	
	TBranch *br_alphaS_MZ = tree->Branch("alphaS_MZ",&alphaS_MZ);
	TBranch *br_alpha = tree->Branch("alpha",&alpha);
	TBranch *br_GFermi = tree->Branch("GFermi",&GFermi);
	TBranch *br_MZ = tree->Branch("MZ",&MZ);
	TBranch *br_Mh = tree->Branch("Mh",&Mh);
	TBranch *br_Mt = tree->Branch("Mt",&Mt);
	TBranch *br_mbmb = tree->Branch("mbmb",&mbmb);
	TBranch *br_Delta_alpha = tree->Branch("Delta_alpha",&Delta_alpha);
	TBranch *br_g = tree->Branch("g",&g);
	TBranch *br_gp = tree->Branch("gp",&gp);
	TBranch *br_g3 = tree->Branch("g3",&g3);
	TBranch *br_yt = tree->Branch("yt",&yt);
	TBranch *br_yb = tree->Branch("yb",&yb);
	TBranch *br_lambda = tree->Branch("lambda",&lambda);
	TBranch *br_x1 = tree->Branch("x1",&x1);
	TBranch *br_x2 = tree->Branch("x2",&x2);
	TBranch *br_x3 = tree->Branch("x3",&x3);
	TBranch *br_y1 = tree->Branch("y1",&y1);
	TBranch *br_y2 = tree->Branch("y2",&y2);
	TBranch *br_z = tree->Branch("z",&z);
	
	int NofsigmasS, NofsigmasZ, Nofsigmash, Nofsigmast, Nperaxis;
	cout << "\nHow many sigmas should be the range for alphaS? ";
	cin >> NofsigmasS;
	cout << "\nHow many sigmas should be the range for MZ? ";
	cin >> NofsigmasZ;
	cout << "\nHow many sigmas should be the range for Mh? ";
	cin >> Nofsigmash;
	cout << "\nHow many sigmas should be the range for Mt? ";
	cin >> Nofsigmast;
	cout << "\nHow many points should be per each axis? ";
	cin >> Nperaxis;
	
	int counter = 0;
	
	for (int iS = 0; iS < Nperaxis; iS++){
	//for (int iGF = 0; iGF < 1; iGF++){
	for (int iZ = 0; iZ < Nperaxis; iZ++){
	for (int ih = 0; ih < Nperaxis; ih++){
	for (int it = 0; it < Nperaxis; it++){
	//for (int ib = 0; ib < 1; ib++){
	//for (int ia = 0; ia < 1; ia++){
		
		alphaS_MZ = SMDR_alphaS_MZ_EXPT + NofsigmasS * SMDR_alphaS_MZ_EXPT_UNC * ( 2 * float(iS) / (Nperaxis-1) - 1);
		alpha = SMDR_alpha_EXPT;
		GFermi = SMDR_GFermi_EXPT;
		MZ = SMDR_MZ_EXPT + NofsigmasZ * SMDR_MZ_EXPT_UNC * ( 2 * float(iZ) / (Nperaxis-1) - 1);
		Mh = SMDR_Mh_EXPT + Nofsigmash * SMDR_Mh_EXPT_UNC * ( 2 * float(ih) / (Nperaxis-1) - 1);
		Mt = SMDR_Mt_EXPT + Nofsigmast * SMDR_Mt_EXPT_UNC * ( 2 * float(it) / (Nperaxis-1) - 1);
		mbmb = SMDR_mbmb_EXPT;
		Delta_alpha = SMDR_Delta_alpha_had_5_MZ_EXPT;	
		
		SMDR_Q_in = 200;
		
		my_Fit_Inputs_original (SMDR_Q_in,
	           alphaS_MZ,
	           alpha,
	           GFermi,
	           MZ,
	           Mh,
	           Mt,
	           mbmb,
	           Delta_alpha,
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
		
		counter += 1;
		cout << "Processed event: " << counter << "/" << Nperaxis*Nperaxis*Nperaxis*Nperaxis/**Nperaxis*Nperaxis*Nperaxis*/ << "\r";
		cout.flush();
				
	}}}}//}}}
	
	tree->Write(("tree_" + to_string(input_config) + "_Npx:" + to_string(Nperaxis) + "_NS:" + to_string(NofsigmasS) + "_NZ:" + to_string(NofsigmasZ) + "_Nh:" + to_string(Nofsigmash) + "_Nt:" + to_string(Nofsigmast)).c_str());
	
	
	file->Close();
	
	}
	
	return 0;

}
