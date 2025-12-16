#include "src/include.h"
#include <stdio.h>

/*-------------------------------------------------------- */
/* Calculation of the abundance of the elements from BBN   */
/*-------------------------------------------------------- */

int main(int argc,char** argv)
{ 
	struct relicparam paramrelic;
	double ratioH[NNUC+1],cov_ratioH[NNUC+1][NNUC+1];
	double H2_H,He3_H,Yp,Li7_H,Li6_H,Be7_H;
	double sigma_H2_H,sigma_He3_H,sigma_Yp,sigma_Li7_H,sigma_Li6_H,sigma_Be7_H;
	double dd0,ndd,Tdrh, dd0_rad;
	double sd0,nsd,Tsrh;
	int coupD;
	int failsafe;
    double Nnu,dNnu,xinu1,xinu2,xinu3;
    double eta,tau,tau_err;
    FILE *fp;
    char DataDir[512];
    char FileNameBase[512];
    char filename[512];

	if(argc<12) 
  	{ 
    printf(" This program needs at least 11 parameters:\n"
		"   dd0     dark energy proportion to photon density at T_i=27e9 K\n"
		"   ndd     dark energy decrease exponent (preferentially >4)\n"
		"   sd0     dark entropy proportion to photon entropy at T=1 MeV\n"
		"   nsd     dark entropy decrease exponent\n"
        "   Nnu    number of neutrinos\n"
        "   dNnu   number of extra neutrino species\n"
        "   dd0_rad  kappa_i of extra radiation species which redshift as a^{-4}\n"
        "   eta    value of the baryon-to-photon ratio\n"
	"   tau      neutron lifetime\n"
        "   DataDir    \n"
        "   FileNameBase    \n");
    printf(" Auxiliary parameters are:\n"
           	"   coup    switch coupling to the plasma (0 or 1, default 1)\n"    
           	"   Td      dark energy cut temperature (in MeV)\n"    
           	"   Ts      dark entropy cut temperature (in MeV)\n"
		"   failsafe    0=fast, 1=precise, 6=robust but slow. See stand_cosmo.c for more options.\n"
                "   tau_err  neutron lifetime uncertainty (optional, needed for error calculation)\n");
      		exit(1); 
  	} 
	else 
	{
        sscanf(argv[1],"%lf",&dd0);
        sscanf(argv[2],"%lf",&ndd);
        sscanf(argv[3],"%lf",&sd0);
        sscanf(argv[4],"%lf",&nsd);
        sscanf(argv[5],"%lf",&Nnu);
        sscanf(argv[6],"%lf",&dNnu);
        sscanf(argv[7],"%lf",&dd0_rad);
        sscanf(argv[8],"%lf",&eta);
	sscanf(argv[9],"%lf",&tau);
        sscanf(argv[10],"%s",DataDir);
        sscanf(argv[11],"%s",FileNameBase);
        if(argc>12) sscanf(argv[12],"%d",&coupD); else coupD=1;		
        if(argc>13) sscanf(argv[13],"%lf",&Tdrh); else Tdrh=0.;
        if(argc>14) sscanf(argv[14],"%lf",&Tsrh); else Tsrh=0.;   
        if(argc>15) sscanf(argv[15],"%d",&failsafe); else failsafe=1;
	if(argc>16) sscanf(argv[16],"%lf",&tau_err); else tau_err=0.5; // Neutron lifetime uncertainty (PDG2022)
  	}
	
	Init_cosmomodel(&paramrelic);	
	
	if(coupD==0) paramrelic.coupd=0; else paramrelic.coupd=1;
	
	paramrelic.failsafe=failsafe;
    
	Init_cosmomodel_param(eta,Nnu,dNnu,dd0_rad,tau,tau_err,0.,0.,0.,&paramrelic);

	Init_dark_density(dd0,ndd,Tdrh*1.e-3,dd0_rad,&paramrelic);
	Init_dark_entropy(sd0,nsd,Tsrh*1.e-3,&paramrelic);



	paramrelic.err=0;
	nucl(&paramrelic,ratioH);
	H2_H=ratioH[3];Yp=ratioH[6];Li7_H=ratioH[8];Be7_H=ratioH[9];He3_H=ratioH[5];Li6_H=ratioH[7];
// 	printf("\t Yp\t\t H2/H\t\t He3/H\t\t Li7/H\t\t Li6/H\t\t Be7/H\n");
//     printf(" \t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\n",Yp,H2_H,He3_H,Li7_H,Li6_H,Be7_H); 


    sprintf(filename, "%s/abundances_%s.txt", DataDir, FileNameBase);

    fp = fopen(filename, "w");
    fprintf(fp, "# Parameters used: \t dd0\t\t ndd\t\t sd0\t\t nsd\t\t Nnu\t\t dNnu\t\t dd0_rad\t\t eta\t\t tau\n");
    fprintf(fp, "# \t\t %.6e\t\t %2.2f\t %.6e\t\t %2.2f\t\t %2.4f\t\t %2.4f\t\t %.6e\t %.6e\t\t %6.2f\n",
	    dd0,ndd,sd0,nsd,Nnu,dNnu,dd0_rad,eta,tau);
    fprintf(fp, "# \t Yp\t\t H2/H\t\t He3/H\t\t Li7/H\t\t Li6/H\t\t Be7/H\n");
    fprintf(fp, " \t %.6e\t %.6e\t %.6e\t %.6e\t %.6e\t %.6e\n",Yp,H2_H,He3_H,Li7_H,Li6_H,Be7_H); 
    fclose(fp);
    


// 	paramrelic.err=3;
// 	if(nucl_err(&paramrelic,ratioH,cov_ratioH))
// 	{
// //         fp = fopen(filename, "w");
// 		fprintf(fp, "# --------------------\n");
// 		fprintf(fp, "# With uncertainties:\n");
//         H2_H=ratioH[3];Yp=ratioH[6];Li7_H=ratioH[8];Be7_H=ratioH[9];He3_H=ratioH[5];Li6_H=ratioH[7];
// 		sigma_H2_H=sqrt(cov_ratioH[3][3]);sigma_Yp=sqrt(cov_ratioH[6][6]);sigma_Li7_H=sqrt(cov_ratioH[8][8]);sigma_Be7_H=sqrt(cov_ratioH[9][9]);sigma_He3_H=sqrt(cov_ratioH[5][5]);sigma_Li6_H=sqrt(cov_ratioH[7][7]);
// 		fprintf(fp, "# \t Yp\t\t H2/H\t\t He3/H\t\t Li7/H\t\t Li6/H\t\t Be7/H\n");
		
// 		fprintf(fp, "# value:\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\n",Yp,H2_H,He3_H,Li7_H,Li6_H,Be7_H); 
// 		fprintf(fp,"# +/- :\n \t %.6e\t %.6e\t %.6e\t %.6e\t %.6e\t %.6e\n\n",sigma_Yp,sigma_H2_H,sigma_He3_H,sigma_Li7_H,sigma_Li6_H,sigma_Be7_H);

// 		double corr_ratioH[NNUC+1][NNUC+1];
// 		for(int ie=1;ie<=NNUC;ie++) for(int je=1;je<=NNUC;je++) corr_ratioH[ie][je]=cov_ratioH[ie][je]/sqrt(cov_ratioH[ie][ie]*cov_ratioH[je][je]);
// 		fprintf(fp, "#Correlation matrix:\n");
// 		fprintf(fp, "#\t Yp\t\t H2/H\t\t He3/H\t\t Li7/H\t\t Li6/H\t\t Be7/H\n");
// 		fprintf(fp, "#Yp\t %f\t %f\t %f\t %f\t %f\t %f\n",corr_ratioH[6][6],corr_ratioH[6][3],corr_ratioH[6][5],corr_ratioH[6][8],corr_ratioH[6][7],corr_ratioH[6][9]);
// 		fprintf(fp, "#H2/H\t %f\t %f\t %f\t %f\t %f\t %f\n",corr_ratioH[3][6],corr_ratioH[3][3],corr_ratioH[3][5],corr_ratioH[3][8],corr_ratioH[3][7],corr_ratioH[3][9]);
// 		fprintf(fp, "#He3/H\t %f\t %f\t %f\t %f\t %f\t %f\n",corr_ratioH[5][6],corr_ratioH[5][3],corr_ratioH[5][5],corr_ratioH[5][8],corr_ratioH[5][7],corr_ratioH[5][9]);
// 		fprintf(fp, "#Li7/H\t %f\t %f\t %f\t %f\t %f\t %f\n",corr_ratioH[8][6],corr_ratioH[8][3],corr_ratioH[8][5],corr_ratioH[8][8],corr_ratioH[8][7],corr_ratioH[8][9]);
// 		fprintf(fp, "#Li6/H\t %f\t %f\t %f\t %f\t %f\t %f\n",corr_ratioH[7][6],corr_ratioH[7][3],corr_ratioH[7][5],corr_ratioH[7][8],corr_ratioH[7][7],corr_ratioH[7][9]);
// 		fprintf(fp, "#Be7/H\t %f\t %f\t %f\t %f\t %f\t %f\n\n",corr_ratioH[9][6],corr_ratioH[9][3],corr_ratioH[9][5],corr_ratioH[9][8],corr_ratioH[9][7],corr_ratioH[9][9]);
// 	}
// 	else fprintf(fp, "#Uncertainty calculation failed\n\n");    

// 	paramrelic.err=2;
// 	nucl(&paramrelic,ratioH);
// 	H2_H=ratioH[3];Yp=ratioH[6];Li7_H=ratioH[8];Be7_H=ratioH[9];He3_H=ratioH[5];Li6_H=ratioH[7];
// 	fprintf(fp, "#  low: \n\t %.6e\t %.6e\t %.6e\t %.6e\t %.6e\t %.6e\n",Yp,H2_H,He3_H,Li7_H,Li6_H,Be7_H);
    
// 	paramrelic.err=1;
// 	nucl(&paramrelic,ratioH);
// 	H2_H=ratioH[3];Yp=ratioH[6];Li7_H=ratioH[8];Be7_H=ratioH[9];He3_H=ratioH[5];Li6_H=ratioH[7];
// 	fprintf(fp, "# high: \n\t %.6e\t %.6e\t %.6e\t %.6e\t %.6e\t %.6e\n\n",Yp,H2_H,He3_H,Li7_H,Li6_H,Be7_H);

    
// 	paramrelic.err=4;
// 	if(nucl_err(&paramrelic,ratioH,cov_ratioH))
// 	{
// 		fprintf(fp, "#--------------------\n\n");
// 		fprintf(fp, "#With MC uncertainties:\n");
//         H2_H=ratioH[3];Yp=ratioH[6];Li7_H=ratioH[8];Be7_H=ratioH[9];He3_H=ratioH[5];Li6_H=ratioH[7];
// 		sigma_H2_H=sqrt(cov_ratioH[3][3]);sigma_Yp=sqrt(cov_ratioH[6][6]);sigma_Li7_H=sqrt(cov_ratioH[8][8]);sigma_Be7_H=sqrt(cov_ratioH[9][9]);sigma_He3_H=sqrt(cov_ratioH[5][5]);sigma_Li6_H=sqrt(cov_ratioH[7][7]);
// 		fprintf(fp, "#\t Yp\t\t H2/H\t\t He3/H\t\t Li7/H\t\t Li6/H\t\t Be7/H\n");
		
// 		fprintf(fp,  "#mean: \n\t %.6e\t %.6e\t %.6e\t %.6e\t %.6e\t %.6e\n",Yp,H2_H,He3_H,Li7_H,Li6_H,Be7_H); 
// 		fprintf(fp, "# +/- : \n\t %.6e\t %.6e\t %.6e\t %.6e\t %.6e\t %.6e\n\n",sigma_Yp,sigma_H2_H,sigma_He3_H,sigma_Li7_H,sigma_Li6_H,sigma_Be7_H);

// 		double corr_ratioH[NNUC+1][NNUC+1];
// 		for(int ie=1;ie<=NNUC;ie++) for(int je=1;je<=NNUC;je++) corr_ratioH[ie][je]=cov_ratioH[ie][je]/sqrt(cov_ratioH[ie][ie]*cov_ratioH[je][je]);
// 		fprintf(fp, "#Correlation matrix:\n");
// 		fprintf(fp, "#\t Yp\t\t H2/H\t\t He3/H\t\t Li7/H\t\t Li6/H\t\t Be7/H\n");
// 		fprintf(fp, "#Yp\t %f\t %f\t %f\t %f\t %f\t %f\n",corr_ratioH[6][6],corr_ratioH[6][3],corr_ratioH[6][5],corr_ratioH[6][8],corr_ratioH[6][7],corr_ratioH[6][9]);
// 		fprintf(fp, "#H2/H\t %f\t %f\t %f\t %f\t %f\t %f\n",corr_ratioH[3][6],corr_ratioH[3][3],corr_ratioH[3][5],corr_ratioH[3][8],corr_ratioH[3][7],corr_ratioH[3][9]);
// 		fprintf(fp, "#He3/H\t %f\t %f\t %f\t %f\t %f\t %f\n",corr_ratioH[5][6],corr_ratioH[5][3],corr_ratioH[5][5],corr_ratioH[5][8],corr_ratioH[5][7],corr_ratioH[5][9]);
// 		fprintf(fp, "#Li7/H\t %f\t %f\t %f\t %f\t %f\t %f\n",corr_ratioH[8][6],corr_ratioH[8][3],corr_ratioH[8][5],corr_ratioH[8][8],corr_ratioH[8][7],corr_ratioH[8][9]);
// 		fprintf(fp, "#Li6/H\t %f\t %f\t %f\t %f\t %f\t %f\n",corr_ratioH[7][6],corr_ratioH[7][3],corr_ratioH[7][5],corr_ratioH[7][8],corr_ratioH[7][7],corr_ratioH[7][9]);
// 		fprintf(fp, "#Be7/H\t %f\t %f\t %f\t %f\t %f\t %f\n\n",corr_ratioH[9][6],corr_ratioH[9][3],corr_ratioH[9][5],corr_ratioH[9][8],corr_ratioH[9][7],corr_ratioH[9][9]);
// 	}
// 	else fprintf(fp, "#Uncertainty calculation failed\n\n");
    
//     fclose(fp);
    
// // 	paramrelic.err=2;
// // 	nucl(&paramrelic,ratioH);
// // 	H2_H=ratioH[3];Yp=ratioH[6];Li7_H=ratioH[8];Be7_H=ratioH[9];He3_H=ratioH[5];Li6_H=ratioH[7];
// // 	printf("  low:\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\n",Yp,H2_H,He3_H,Li7_H,Li6_H,Be7_H);
    
// //	printf(" cent:\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\n",Yp,H2_H,He3_H,Li7_H,Li6_H,Be7_H); 
    
    
// // 	paramrelic.err=1;
// // 	nucl(&paramrelic,ratioH);
// // 	H2_H=ratioH[3];Yp=ratioH[6];Li7_H=ratioH[8];Be7_H=ratioH[9];He3_H=ratioH[5];Li6_H=ratioH[7];
// // 	printf(" high:\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\n\n",Yp,H2_H,He3_H,Li7_H,Li6_H,Be7_H);
			
// // 	paramrelic.err=3;
// // 	if(nucl_err(&paramrelic,ratioH,cov_ratioH))
// // 	{
// // 		printf("--------------------\n\n");
// // 		printf("With uncertainties:\n");
// //         H2_H=ratioH[3];Yp=ratioH[6];Li7_H=ratioH[8];Be7_H=ratioH[9];He3_H=ratioH[5];Li6_H=ratioH[7];
// // 		sigma_H2_H=sqrt(cov_ratioH[3][3]);sigma_Yp=sqrt(cov_ratioH[6][6]);sigma_Li7_H=sqrt(cov_ratioH[8][8]);sigma_Be7_H=sqrt(cov_ratioH[9][9]);sigma_He3_H=sqrt(cov_ratioH[5][5]);sigma_Li6_H=sqrt(cov_ratioH[7][7]);
// // 		printf("\t Yp\t\t H2/H\t\t He3/H\t\t Li7/H\t\t Li6/H\t\t Be7/H\n");
		
// // 		printf("value:\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\n",Yp,H2_H,He3_H,Li7_H,Li6_H,Be7_H); 
// // 		printf(" +/- :\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\n\n",sigma_Yp,sigma_H2_H,sigma_He3_H,sigma_Li7_H,sigma_Li6_H,sigma_Be7_H);

// // 		double corr_ratioH[NNUC+1][NNUC+1];
// // 		for(int ie=1;ie<=NNUC;ie++) for(int je=1;je<=NNUC;je++) corr_ratioH[ie][je]=cov_ratioH[ie][je]/sqrt(cov_ratioH[ie][ie]*cov_ratioH[je][je]);
// // 		printf("Correlation matrix:\n");
// // 		printf("\t Yp\t\t H2/H\t\t He3/H\t\t Li7/H\t\t Li6/H\t\t Be7/H\n");
// // 		printf("Yp\t %f\t %f\t %f\t %f\t %f\t %f\n",corr_ratioH[6][6],corr_ratioH[6][3],corr_ratioH[6][5],corr_ratioH[6][8],corr_ratioH[6][7],corr_ratioH[6][9]);
// // 		printf("H2/H\t %f\t %f\t %f\t %f\t %f\t %f\n",corr_ratioH[3][6],corr_ratioH[3][3],corr_ratioH[3][5],corr_ratioH[3][8],corr_ratioH[3][7],corr_ratioH[3][9]);
// // 		printf("He3/H\t %f\t %f\t %f\t %f\t %f\t %f\n",corr_ratioH[5][6],corr_ratioH[5][3],corr_ratioH[5][5],corr_ratioH[5][8],corr_ratioH[5][7],corr_ratioH[5][9]);
// // 		printf("Li7/H\t %f\t %f\t %f\t %f\t %f\t %f\n",corr_ratioH[8][6],corr_ratioH[8][3],corr_ratioH[8][5],corr_ratioH[8][8],corr_ratioH[8][7],corr_ratioH[8][9]);
// // 		printf("Li6/H\t %f\t %f\t %f\t %f\t %f\t %f\n",corr_ratioH[7][6],corr_ratioH[7][3],corr_ratioH[7][5],corr_ratioH[7][8],corr_ratioH[7][7],corr_ratioH[7][9]);
// // 		printf("Be7/H\t %f\t %f\t %f\t %f\t %f\t %f\n\n",corr_ratioH[9][6],corr_ratioH[9][3],corr_ratioH[9][5],corr_ratioH[9][8],corr_ratioH[9][7],corr_ratioH[9][9]);
// // 	}
// // 	else printf("Uncertainty calculation failed\n\n");

// // 	/*paramrelic.err=4;
// // 	if(nucl_err(&paramrelic,ratioH,cov_ratioH))
// // 	{
// // 		printf("--------------------\n\n");
// // 		printf("With MC uncertainties:\n");
// //         H2_H=ratioH[3];Yp=ratioH[6];Li7_H=ratioH[8];Be7_H=ratioH[9];He3_H=ratioH[5];Li6_H=ratioH[7];
// // 		sigma_H2_H=sqrt(cov_ratioH[3][3]);sigma_Yp=sqrt(cov_ratioH[6][6]);sigma_Li7_H=sqrt(cov_ratioH[8][8]);sigma_Be7_H=sqrt(cov_ratioH[9][9]);sigma_He3_H=sqrt(cov_ratioH[5][5]);sigma_Li6_H=sqrt(cov_ratioH[7][7]);
// // 		printf("\t Yp\t\t H2/H\t\t He3/H\t\t Li7/H\t\t Li6/H\t\t Be7/H\n");
		
// // 		printf("mean:\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\n",Yp,H2_H,He3_H,Li7_H,Li6_H,Be7_H); 
// // 		printf(" +/- :\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\n\n",sigma_Yp,sigma_H2_H,sigma_He3_H,sigma_Li7_H,sigma_Li6_H,sigma_Be7_H);

// // 		double corr_ratioH[NNUC+1][NNUC+1];
// // 		for(int ie=1;ie<=NNUC;ie++) for(int je=1;je<=NNUC;je++) corr_ratioH[ie][je]=cov_ratioH[ie][je]/sqrt(cov_ratioH[ie][ie]*cov_ratioH[je][je]);
// // 		printf("Correlation matrix:\n");
// // 		printf("\t Yp\t\t H2/H\t\t He3/H\t\t Li7/H\t\t Li6/H\t\t Be7/H\n");
// // 		printf("Yp\t %f\t %f\t %f\t %f\t %f\t %f\n",corr_ratioH[6][6],corr_ratioH[6][3],corr_ratioH[6][5],corr_ratioH[6][8],corr_ratioH[6][7],corr_ratioH[6][9]);
// // 		printf("H2/H\t %f\t %f\t %f\t %f\t %f\t %f\n",corr_ratioH[3][6],corr_ratioH[3][3],corr_ratioH[3][5],corr_ratioH[3][8],corr_ratioH[3][7],corr_ratioH[3][9]);
// // 		printf("He3/H\t %f\t %f\t %f\t %f\t %f\t %f\n",corr_ratioH[5][6],corr_ratioH[5][3],corr_ratioH[5][5],corr_ratioH[5][8],corr_ratioH[5][7],corr_ratioH[5][9]);
// // 		printf("Li7/H\t %f\t %f\t %f\t %f\t %f\t %f\n",corr_ratioH[8][6],corr_ratioH[8][3],corr_ratioH[8][5],corr_ratioH[8][8],corr_ratioH[8][7],corr_ratioH[8][9]);
// // 		printf("Li6/H\t %f\t %f\t %f\t %f\t %f\t %f\n",corr_ratioH[7][6],corr_ratioH[7][3],corr_ratioH[7][5],corr_ratioH[7][8],corr_ratioH[7][7],corr_ratioH[7][9]);
// // 		printf("Be7/H\t %f\t %f\t %f\t %f\t %f\t %f\n\n",corr_ratioH[9][6],corr_ratioH[9][3],corr_ratioH[9][5],corr_ratioH[9][8],corr_ratioH[9][7],corr_ratioH[9][9]);
// // 	}
// // 	else printf("Uncertainty calculation failed\n\n");*/
		
// // 	paramrelic.err=0;
// // 	int compat=bbn_excluded(&paramrelic);

// // 	if(compat==1) printf("Excluded by BBN constraints (chi2 without correlations)\n");
// // 	else if(compat==0) printf("Compatible with BBN constraints (chi2 without correlations)\n");
// // 	else printf("Computation failed (chi2 without correlations)\n");

// // 	paramrelic.err=3;
// // 	compat=bbn_excluded(&paramrelic);

// // 	if(compat==1) printf("Excluded by BBN constraints (chi2 including correlations)\n");
// // 	else if(compat==0) printf("Compatible with BBN constraints (chi2 including correlations)\n");
// // 	else printf("Computation failed (chi2 including correlations)\n");

	return 1;
}
