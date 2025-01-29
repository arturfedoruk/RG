// to launch the program: 
// g++ `root-config --cflags` withSMDRgeneraterandom.cpp `root-config --libs` -lm -lsmdr -ltsil -l3vil

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

	

	return 0;

}
