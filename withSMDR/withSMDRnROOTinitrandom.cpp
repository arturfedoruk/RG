// to launch the program: 
// g++ `root-config --cflags` withSMDRnROOTinitrandom.cpp `root-config --libs` -lm -lsmdr -ltsil -l3vil

#include "smdr.h"
#include "iostream"
#include "fstream"
#include "string"
#include "TROOT.h"
#include "TError.h"
#include "TTree.h"
#include "TFile.h"
using namespace std;
#define ZEROSAFE(a) (((a) > (SMDR_TOL)) ? (a) : (SMDR_TOL)) //idk wht's that

#include "my_Fit_Inputs.cpp"
#include "loop_configs.cpp"

int main(){
	
	#include "smdr_pdg_2025.h"

	TFile *file = new TFile("random_data.root","recreate");
	
	const int nconfigs = 6;
	
	TTree *tree_111111 = new TTree("tree_111111","tree_111111");
	TTree *tree_222222 = new TTree("tree_222222","tree_222222");
	TTree *tree_333333 = new TTree("tree_333333","tree_333333");
	TTree *tree_333221 = new TTree("tree_333221","tree_333221");
	TTree *tree_444332 = new TTree("tree_444332","tree_444332");
	TTree *tree_444333 = new TTree("tree_444333","tree_444333");
	
	TTree* trees[nconfigs] = {tree_111111, tree_222222, tree_333333, tree_333221, tree_444332, tree_444333};
	
	char* tree_names[nconfigs] = {"tree_111111", "tree_222222", "tree_333333", "tree_333221", "tree_444332", "tree_444333"};
	
	int seed;
	double alphaS_MZ, alpha, GFermi, MZ, Mh, Mt, mbmb, Delta_alpha, g, gp, g3, yt, yb, lambda, x1, x2, x3, y1, y2, z;
	
	double ERROR_TOLERANCE = 1e-6;
	double pi = 3.14159265359;
	
	for (int itree = 0; itree < nconfigs; itree++){
		SMDR_Q_in = 173.22;
		
		cout << "Creating " << tree_names[itree] << endl;
		
		TBranch *br_seed = trees[itree]->Branch("seed",&seed);
		TBranch *br_alphaS_MZ = trees[itree]->Branch("alphaS_MZ",&alphaS_MZ);
		TBranch *br_alpha = trees[itree]->Branch("alpha",&alpha);
		TBranch *br_GFermi = trees[itree]->Branch("GFermi",&GFermi);
		TBranch *br_MZ = trees[itree]->Branch("MZ",&MZ);
		TBranch *br_Mh = trees[itree]->Branch("Mh",&Mh);
		TBranch *br_Mt = trees[itree]->Branch("Mt",&Mt);
		TBranch *br_mbmb = trees[itree]->Branch("mbmb",&mbmb);
		TBranch *br_Delta_alpha = trees[itree]->Branch("Delta_alpha",&Delta_alpha);
		TBranch *br_g = trees[itree]->Branch("g",&g);
		TBranch *br_gp = trees[itree]->Branch("gp",&gp);
		TBranch *br_g3 = trees[itree]->Branch("g3",&g3);
		TBranch *br_yt = trees[itree]->Branch("yt",&yt);
		TBranch *br_yb = trees[itree]->Branch("yb",&yb);
		TBranch *br_lambda = trees[itree]->Branch("lambda",&lambda);
		TBranch *br_x1 = trees[itree]->Branch("x1",&x1);
		TBranch *br_x2 = trees[itree]->Branch("x2",&x2);
		TBranch *br_x3 = trees[itree]->Branch("x3",&x3);
		TBranch *br_y1 = trees[itree]->Branch("y1",&y1);
		TBranch *br_y2 = trees[itree]->Branch("y2",&y2);
		TBranch *br_z = trees[itree]->Branch("z",&z);
		
		seed = -1;
		alphaS_MZ = SMDR_alphaS_MZ_EXPT;
		alpha = SMDR_alpha_EXPT;
		GFermi = SMDR_GFermi_EXPT;
		MZ = SMDR_MZ_EXPT;
		Mh = SMDR_Mh_EXPT;
		Mt = SMDR_Mt_EXPT;
		mbmb = SMDR_mbmb_EXPT;
		Delta_alpha = SMDR_Delta_alpha_had_5_MZ_EXPT;
		
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
		           loop_configs[itree]);
		           
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
		
		trees[itree]->Fill();
		trees[itree]->Write(tree_names[itree]);
	
	}
	
	file->Close();
	
	ofstream lastSeed("lastSeed.txt");
	lastSeed << 0 << endl;
	lastSeed.close();
	
	return 0;

}
