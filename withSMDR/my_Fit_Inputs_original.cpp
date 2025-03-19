// version of SMDR_Fit_inputs with added loop customization and without particles lighter than b (all Qs are original)
int  my_Fit_Inputs_original (SMDR_REAL Q_target,
                      SMDR_REAL alphaS_5_MZ_target,
                      SMDR_REAL alpha_target,
                      SMDR_REAL GFermi_target,
                      SMDR_REAL MZ_target,
                      SMDR_REAL Mh_target,
                      SMDR_REAL Mt_target,
                      SMDR_REAL mbmb_target,
                       // i totally eliminated all leptons and quarks lighter than b from the whole code (they're too smol to account for)
                      SMDR_REAL Delta_alpha_had_target,
                      SMDR_REAL error_target,
                      float* loop_config)
{
  char funcname[] = "my_Fit_Inputs";
  int j;
  SMDR_REAL vratio, v2ratio;
  SMDR_REAL KZ, Ka;
  SMDR_REAL err_total, err_alphaS, err_alpha, err_MZ, err_GFermi, err_Mh; 
  SMDR_REAL err_Mt, err_mbmb;  
  SMDR_REAL g3rat, lambdarat, ytrat, ybrat;
  SMDR_REAL del_v, del_g3, del_g, del_gp, del_yt, del_yb, del_lambda;
  int result_niters;
  
  // accelcoefs.h literally
	SMDR_REAL cog3_g3     = 0.9098046982;
	SMDR_REAL cog3_g32    = -0.1287280724;
	SMDR_REAL colambda_g3 = -0.0177478612;
	SMDR_REAL coyt_g3     = -0.1604767897;
	SMDR_REAL coyb_g3     = -1.4310376608;
	SMDR_REAL coyt_yt     = 1.0872686457;
	SMDR_REAL coyt_yt2    = 0.0127850017;
	SMDR_REAL colambda_yt = 0.1002393382;
	SMDR_REAL cog3_yt     = -0.0069382506;
	SMDR_REAL coyb_yt     = -0.0070378697;
	SMDR_REAL coyb_yb     = 1.1850813952;
	SMDR_REAL coyb_yb2    = 0.0744955866;

  if (mbmb_target < 3.0) {
    SMDR_Warn (funcname, "Target mb(mb) is too small.");
    SMDR_Warn (funcname, "Setting it and all lighter fermion masses to 0.\n");
    mbmb_target = 0;
  } 

  /* Don't compute Mt, Mh with needlessly high tolerances... */ 
  SMDR_MTPOLE_TOLERANCE = 1.0e-8;
  SMDR_MHPOLE_TOLERANCE = 1.0e-8;
  if (error_target < 1.0e-7) {
    SMDR_MTPOLE_TOLERANCE = error_target/10.;
    SMDR_MHPOLE_TOLERANCE = error_target/10.;
  }
  if (SMDR_MTPOLE_TOLERANCE < 1.01e-15) SMDR_MTPOLE_TOLERANCE = 1.01e-15;
  if (SMDR_MHPOLE_TOLERANCE < 1.01e-15) SMDR_MHPOLE_TOLERANCE = 1.01e-15;

  /* Use ReferenceModel.dat as our initial guess for the MSbar
     parameters, and their corresponding OS predictions.
  */
  /* SMDR_Read_MSbar_Inputs ("ReferenceModel.dat"); */
  /* SMDR_Read_Value ("ReferenceModel.dat", "SMDR_v_in"); */
  /* SMDR_Read_OS_Inputs ("ReferenceModel.dat"); */
  
	// AUTOrefmodel.h literally	
	SMDR_Q_in  = 172.5;  
	SMDR_g3_in = 1.1630062487549031;
	SMDR_g_in  = 0.6476067573057625;
	SMDR_gp_in = 0.3585502116947177;
	SMDR_v_in  = 246.60321691251776;
	SMDR_lambda_in  = 0.1263927658461975;
	SMDR_yt_in = 0.9315770153451640;
	SMDR_yb_in = 0.0155032393874445;
	SMDR_Delta_alpha_had_5_MZ_in = 0.02766;
	SMDR_m2_in = -8633.496280611358;
	SMDR_Lambda_in = 123709045.022956;
	SMDR_Mt_pole = 172.5;
	SMDR_Mh_pole = 125.25;
	SMDR_MZ_PDG = 91.1876;
	SMDR_MZ_pole = 91.16195836913282;
	SMDR_MW_PDG = 80.35247573675088;
	SMDR_MW_pole = 80.33209076550887;
	SMDR_GFermi = 0.000011663787;
	SMDR_alpha = 1/137.0359990836963270;
	SMDR_alphaS_5_MZ = 0.1179;
	SMDR_mbmb = 4.18;

  SMDR_RGeval_SM (Mt_target, loop_config[6]);

  SMDR_Save_Inputs ();
  SMDR_Delta_alpha_had_5_MZ_in = Delta_alpha_had_target;

  for (j=0; j<18; j++) {

    /* Find the new input parameter estimates by comparing OS observables to the 
       targets. On the 0th iteration, the OS observables and MSbar parameters
       are both as specified in "ReferenceModel.dat".
    */
    v2ratio = SMDR_GFermi/GFermi_target;
    vratio = SMDR_SQRT (v2ratio);

    KZ = (SMDR_g_in * SMDR_g_in + SMDR_gp_in * SMDR_gp_in) * 
         (MZ_target * MZ_target)/
         (SMDR_MZ_PDG * SMDR_MZ_PDG * v2ratio);

    Ka = ((SMDR_g_in * SMDR_g_in * SMDR_gp_in * SMDR_gp_in)/
          (SMDR_g_in * SMDR_g_in + SMDR_gp_in * SMDR_gp_in)) * 
          (alpha_target/SMDR_alpha);
	
	//here I made g_in -> SMDR_g_in and gp_in -> SMDR_gp_in
    del_g  = (SMDR_SQRT(KZ * (1. + SMDR_SQRT(1. - 4. * Ka/KZ))/2.)/SMDR_g_in) - 1.;
    del_gp = (SMDR_SQRT(KZ * (1. - SMDR_SQRT(1. - 4. * Ka/KZ))/2.)/SMDR_gp_in) - 1.;

    del_g3 = SMDR_SQRT (alphaS_5_MZ_target/SMDR_alphaS_5_MZ) - 1.;
    del_v = vratio - 1.;
    del_lambda = ((Mh_target * Mh_target)/
                 (SMDR_Mh_pole * SMDR_Mh_pole * v2ratio)) - 1.;
    del_yt = (Mt_target/(SMDR_Mt_pole * vratio)) - 1.;
    del_yb = (mbmb_target/(ZEROSAFE(SMDR_mbmb) * vratio)) - 1.;

    g3rat = 1. + cog3_g3 * del_g3 
               + cog3_g32 * del_g3 * del_g3 
               + cog3_yt * del_yt;

    lambdarat = 1. + del_lambda 
                   + colambda_yt * del_yt 
                   + colambda_g3 * del_g3;

    ytrat = 1. + coyt_yt * del_yt 
               + coyt_yt2 * del_yt * del_yt
               + coyt_g3 * del_g3;

    ybrat = 1. + coyb_yb * del_yb 
               + coyb_yb2 * del_yb * del_yb +
               + coyb_g3 * del_g3
               + coyb_yt * del_yt;

    if ((ybrat < 0) || (mbmb_target < SMDR_TOL)) ybrat = 0;

    SMDR_v_in *= vratio;
    SMDR_g3_in *= g3rat;
    SMDR_g_in  *= 1. + del_g;
    SMDR_gp_in *= 1. + del_gp;
    SMDR_lambda_in *= lambdarat;
    SMDR_yt_in *= ytrat;
    SMDR_yb_in *= ybrat;

    err_alphaS = SMDR_alphaS_5_MZ/alphaS_5_MZ_target - 1.;
    err_alpha = SMDR_alpha/alpha_target - 1.;
    err_MZ = SMDR_MZ_PDG/MZ_target - 1.;
    err_GFermi = SMDR_GFermi/GFermi_target - 1.;
    err_Mh = SMDR_Mh_pole/Mh_target - 1.;
    err_Mt = SMDR_Mt_pole/Mt_target - 1.;
    err_mbmb = (SMDR_mbmb - mbmb_target)/SMDR_mbmb_EXPT;

    err_total = SMDR_FABS(err_alphaS) +
                SMDR_FABS(err_alpha) + 
                SMDR_FABS(err_MZ) + 
                SMDR_FABS(err_GFermi) + 
                SMDR_FABS(err_Mh) + 
                SMDR_FABS(err_Mt) + 
                SMDR_FABS(err_mbmb);
                
      // i deleted some commented sectors for debugging, for full version see fit_obc.c

    if (j>0) {
      /*printf("Iteration %d: ",j);
      printf("Total fractional error = %.15Lf\n", err_total); */
      if (err_total < error_target) break;
    }
  
    /* Now evaluate all of the on-shell observables. */
    SMDR_Eval_Mt_pole (SMDR_Mt_EXPT, 1, loop_config[0], loop_config[1], &SMDR_Mt_pole, &SMDR_Gammat_pole);
    SMDR_Eval_Mh_pole (160., loop_config[2], &SMDR_Mh_pole, &SMDR_Gammah_pole);
    SMDR_Eval_MZ_pole (160., loop_config[3], &SMDR_MZ_pole, &SMDR_GammaZ_pole,
                       &SMDR_MZ_PDG, &SMDR_GammaZ_PDG);

    /* This one is computed only because it is needed by SMDR_Eval_Gauge() */
    SMDR_Eval_MW_pole (160., loop_config[4], &SMDR_MW_pole, &SMDR_GammaW_pole,
                       &SMDR_MW_PDG, &SMDR_GammaW_PDG);

    SMDR_GFermi = SMDR_Eval_GFermi (SMDR_Mt_EXPT, loop_config[5]);
    SMDR_Eval_Gauge (SMDR_Mt_pole, SMDR_Mh_pole, SMDR_MW_PDG);
    SMDR_Eval_QCDQED_at_MZ (SMDR_MZ_EXPT, SMDR_MZ_EXPT, loop_config[6]);

    if (mbmb_target > SMDR_TOL) 
      SMDR_mbmb = SMDR_Eval_mbmb(SMDR_MZ_EXPT, loop_config[7]);
    else SMDR_mbmb = 0;

    result_niters = j;
  }

  SMDR_Load_Inputs();

  SMDR_Lambda = 0;
  SMDR_Lambda_in = SMDR_Lambda = -SMDR_Eval_Veffmin (-1, loop_config[2]);
  SMDR_m2_in = SMDR_m2;

  SMDR_RGeval_SM (Q_target, loop_config[6]);
  SMDR_Save_Inputs();

  return(result_niters);
};
