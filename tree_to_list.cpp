#include "iostream"
#include "fstream"
#include "string"
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
using namespace std;

TFile *file;
TTree *tree;

int main(){

	string filename_temp, treename_temp;
	cout << "Enter file name " << endl;
	cin >> filename_temp;
	cout << "Enter name of tree (note that it must have on-shell and MSbar branches!) " << endl;
	cin >> treename_temp;
	
	const char* filename = filename_temp.c_str();
	const char* treename = treename_temp.c_str();
	
	file = new TFile(filename);
	tree = (TTree*)file->Get(treename);
	
	if (!tree) cout << "TREE NOT FOUND" << endl;
	
	ofstream outfile(("converted_" + treename_temp + ".txt").c_str(), ios::app | ios::out);
	outfile << "{" << endl << "{ {alphaS, Mt, Mh, MZ, Delta, mbmb, GF}, {g, gp, g3, yt, yb, lambda} }," << endl;
	
	double alphaS, Mt, Mh, MZ, Delta, mbmb, GF, g, gp, g3, yt, yb, lambda;
	
	tree->SetBranchAddress("g",&g);
	tree->SetBranchAddress("gp",&gp);
	tree->SetBranchAddress("g3",&g3);
	tree->SetBranchAddress("yt",&yt);
	tree->SetBranchAddress("yb",&yb);
	tree->SetBranchAddress("lambda",&lambda);
	tree->SetBranchAddress("alphaS_MZ",&alphaS);
	tree->SetBranchAddress("Mt",&Mt);
	tree->SetBranchAddress("Mh",&Mh);
	tree->SetBranchAddress("MZ",&MZ);
	tree->SetBranchAddress("Delta_alpha",&Delta);
	tree->SetBranchAddress("mbmb",&mbmb);
	tree->SetBranchAddress("GFermi",&GF);
	
	int N = tree->GetEntries();
	for(int ientry = 0; ientry < N; ientry++){
			
		tree->GetEntry(ientry);
		
		outfile << "{ {" << alphaS << ", " << Mt << ", " << Mh << ", " << MZ << ", " << Delta << ", " << mbmb << ", " << GF << "}, {" << g << ", " << gp << ", " << g3 << ", " << yt << ", " << yb << ", " << lambda << "} }," << endl;

		// counter (idk, the process is quite fast, but i added this)
		cout << "Processed event: " << ientry+1 << "/" << N << "\r";
		cout.flush();
		
		}
	
	outfile << "}," << endl;
	return 0;

}
