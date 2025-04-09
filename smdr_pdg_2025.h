// data to be used in paper, most recent values 
// PDG_ values are used in case SMDR is not included (e.g. fit.cpp)

double PDG_GFermi =     0.000011663787;
double PDG_GFermi_UNC = 0.000000000006;

double PDG_alpha =     0.0072973525693; 
double PDG_alpha_UNC = 0.0000000000011;

double PDG_alphaS_MZ =     0.1180;     
double PDG_alphaS_MZ_UNC = 0.0009;

double PDG_Delta_alpha_had_5_MZ =     0.02783;
double PDG_Delta_alpha_had_5_MZ_UNC = 0.00006;

double PDG_Mt =   172.57;    
double PDG_Mt_UNC = 0.29;

double PDG_Mh =   125.20;
double PDG_Mh_UNC = 0.11;

double PDG_MZ =    91.1880; 
double PDG_MZ_UNC = 0.0020;

double PDG_MW =    80.3692;  
double PDG_MW_UNC = 0.0133;

double PDG_mbmb =        4.183;   
double PDG_mbmb_UNC_hi = 0.007;

#ifdef _SMDR_H_

SMDR_GFermi_EXPT =     PDG_GFermi;
SMDR_GFermi_EXPT_UNC = PDG_GFermi_UNC;

SMDR_alpha_EXPT =     PDG_alpha; 
SMDR_alpha_EXPT_UNC = PDG_alpha_UNC;

SMDR_alphaS_MZ_EXPT =     PDG_alphaS_MZ;     
SMDR_alphaS_MZ_EXPT_UNC = PDG_alphaS_MZ_UNC;

SMDR_Delta_alpha_had_5_MZ_EXPT =     PDG_Delta_alpha_had_5_MZ;
SMDR_Delta_alpha_had_5_MZ_EXPT_UNC = PDG_Delta_alpha_had_5_MZ_UNC;

SMDR_Mt_EXPT =   PDG_Mt;    
SMDR_Mt_EXPT_UNC = PDG_Mt_UNC;

SMDR_Mh_EXPT =   PDG_Mh;
SMDR_Mh_EXPT_UNC = PDG_Mh_UNC;

SMDR_MZ_EXPT =    PDG_MZ; 
SMDR_MZ_EXPT_UNC = PDG_MZ_UNC;

SMDR_MW_EXPT =    PDG_MW;  
SMDR_MW_EXPT_UNC = PDG_MW_UNC;

SMDR_mbmb_EXPT =        PDG_mbmb;   
SMDR_mbmb_EXPT_UNC_hi = PDG_mbmb_UNC_hi;

#endif
