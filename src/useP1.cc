//useP1.cc
#include <cstdlib>
#include <R.h>
#include <new>
#include "p1.h"

using namespace std;
using namespace TNT;

extern "C" {void FGweight_KM(float *ts, int *ics, float *zs, int *N, int *NCOV, float *gfg1, float *var1, float *tbase_km, float *lambda10km, float *lambda10sd_km, float *Wlambda, float *Wbeta, float *conv, float *variance, float *var_conservative);}
extern "C" {void FGweight_COX(float *ts, int *ics, float *zs, float *zcs, int *n, int *ncov, int *ncovc, float *gfg3, float *var3, float *tbase_cox, float *lambda10cox, float *lambda10sd_cox, float *Wlambda, float *Wbeta, float *conv, float *censdet, float *variance, float *var_conservative);}
extern "C" {void FGweight_KM_Strata(float *ts, int *ics, float *zs, int *N, int *NCOV, float *gfg1, float *var1, float *tbase_km, float *lambda10km, float *lambda10sd_km, float *Wlambda, float *Wbeta, float *conv, int *NSTRATA, int *strata_var, float *variance, float *var_conservative);}
extern "C" {void FGweight_COX_Strata(float *ts, int *ics, float *zs, float *zcs, int *N, int *NCOV, int *NCOVC, float *gfg3, float *var3, float *tbase_cox, float *lambda10cox, float *lambda10sd_cox, float *Wlambda, float *Wbeta, float *conv, float *censdet,int *NSTRATA, int *strata_var, float *variance, float *var_conservative);}
				 
// .C interface for "KM"				  
void FGweight_KM(float *ts, int *ics, float *zs, int *N, int *NCOV, float *gfg1, float *var1, float *tbase_km, float *lambda10km, float *lambda10sd_km, float *Wlambda, float *Wbeta, float *conv, float *variance, float *var_conservative) {
	//variable declaration
	int n = *N;
	int ncov = *NCOV;
	int nfgdiv = 0, nit = 0;
	int njp = 0.;
	int *njpp = &njp;
	
	float iflag = 0.;
	float *ifl = &iflag;
	float max_error = 1.;
	
	vector<float> dnc(n);
	vector<float> yc(n);	
	vector<float> tjp(n);
	vector<float> g1(ncov);
	vector<float> g2(ncov);
	vector<float> error(ncov);
	
	float **wdn1 = new float*[n];
	for(int i=0; i<n; i++) wdn1[i] = new float[n];
	float **wy1 = new float*[n];
	for(int i=0; i<n; i++) wy1[i] = new float[n];	
	float **dnci = new float*[n];
	for(int i=0; i<n; i++) dnci[i] = new float[n];
	
	FG_KM W(n,ncov,ts,ics,zs);
	W.riskst(wdn1, wy1, dnc, yc, dnci, njpp, tjp);
	for(int i=0; i<ncov; i++) {
		g1[i] = 0.;
		g2[i] = 0.;
	}
	nit = 0;
	
	//newton-raphson for FG model
	while(max_error > 0.00001 && nit < 50) {
		W.fgest1(njp, wy1, wdn1, g1, g2, ifl, ncov);
		nit++;
		if(iflag==1) {
			nfgdiv++;
			Rprintf("Algorithm did not converge when fitting the model with a Kaplan-Meier weight\n Program stopped.\n\n");
			conv[0] = 0.;
			return;
		}
		else {
			max_error = fabs(g2[0]-g1[0]);
			for(int k=1; k<ncov; k++) {
				error[k] = fabs(g2[k]-g1[k]);		
				if(error[k] > max_error) max_error = error[k];
			}
			for(int k=0; k<ncov; k++) g1[k] = g2[k];		
		}//endif	
	}//end while
	if (nit==50) {
		nfgdiv++;
		Rprintf("Algorithm did not converge when fitting the model with a Kaplan-Meier weight\n Program stopped.\n\n");
		conv[0] = 0.;
		return;
	}	
	for(int k=0; k<ncov; k++) gfg1[k]=g2[k];	
	
	//variance calculation
	if(variance[0] == 1.) {
		W.fgest2km(tjp,njp,wy1,wdn1,gfg1,dnc,yc,dnci,var1,ncov,lambda10km,tbase_km,lambda10sd_km,Wlambda,Wbeta,var_conservative);
	}
	
	//memory cleanup
	for(int i=0; i<n; i++) delete[] wdn1[i];
	delete[] wdn1;
	for(int i=0; i<n; i++) delete[] wy1[i];
	delete[] wy1;
	for(int i=0; i<n; i++) delete[] dnci[i];
	delete[] dnci;
}


//.C interface for "COX"
void FGweight_COX(float *ts, int *ics, float *zs, float *zcs, int *N, int *NCOV, int *NCOVC, float *gfg3, float *var3, float *tbase_cox, float *lambda10cox, float *lambda10sd_cox, float *Wlambda, float *Wbeta, float *conv, float *censdet, float *variance, float *var_conservative) {
	//variable declaration
	int n = *N;
	int ncov = *NCOV;
	int ncovc = *NCOVC;
	int ncsdiv = 0, nfgdiv = 0, nit = 0;
	int njp = 0.;
	int *njpp = &njp;		

	float iflagc = 0.;
	float *iflc = &iflagc;
	float iflag = 0.;
	float *ifl = &iflag;
	float iflagric = 0.;
	float *iflric = &iflagric;
	float max_error = 1.;
	
	vector<float> dnc(n);
	vector<float> yc(n);
	vector<float> tjp(n);
	vector<float> rlamc0(n);
	vector<float> g1(ncov);
	vector<float> g2(ncov);
	vector<float> betac(ncovc);
	vector<float> error(ncov);
	
	float **wdn3 = new float*[n];
	for(int i=0; i<n; i++) wdn3[i] = new float[n];	
	float **wy3 = new float*[n];
	for(int i=0; i<n; i++) wy3[i] = new float[n];		
	float **dnci = new float*[n];
	for(int i=0; i<n; i++) dnci[i] = new float[n];
	float **s1s0c = new float*[n];
	for(int i=0; i<n; i++) s1s0c[i] = new float[n];
	float **ric = new float*[ncovc];
	for(int i=0; i<ncovc; i++) ric[i] = new float[ncovc];
	
	FG_COX W(n,ncov,ncovc,ts,ics,zs,zcs);
	W.riskst(wdn3, wy3, iflc, betac, dnc, yc, dnci, njpp, tjp, ric, rlamc0, s1s0c, ncovc, censdet);
	if(iflagc==1) {
		ncsdiv++;
		Rprintf("Warning message:\n Algorithm did not converge when fitting the censoring distribution.\n Program stopped.\n\n");
		conv[0] = 0.;
		return;
	}
	for(int i=0; i<ncov; i++) {
		g1[i] = 0.;
		g2[i] = 0.;
	}
	nit = 0;

	//newton-raphson for FG model
	while(max_error > 0.00001 && nit < 50) {
		W.fgest1(njp, wy3, wdn3, g1, g2, ifl, ncov);
		nit++;
		if(iflag==1) {
			nfgdiv++;
			Rprintf("Algorithm did not converge when fitting the model with a Cox weight\n Program stopped.\n\n");
			conv[0] = 0.;
			return;
		}
		else {
			max_error = fabs(g2[0]-g1[0]);
			for(int k=1; k<ncov; k++) {
				error[k] = fabs(g2[k]-g1[k]);		
				if(error[k] > max_error) max_error = error[k];
			}
			for(int k=0; k<ncov; k++) g1[k] = g2[k];		
		}//endif	
	}//end while
	if (nit==50) {
		nfgdiv++;
		Rprintf("Algorithm did not converge when fitting the model with a Cox weight\n Program stopped.\n\n");
		conv[0] = 0.;
		return;
	}	
	for(int k=0; k<ncov; k++) gfg3[k]=g2[k];	

	//variance calculation
	if(variance[0] == 1.) {
		W.fgest2cox(tjp,njp,betac,ric,rlamc0,s1s0c,wy3,wdn3,gfg3,dnc,yc,dnci,var3,lambda10cox,tbase_cox,lambda10sd_cox,Wlambda,Wbeta,iflric,var_conservative);
		if (iflagric==1) {
			Rprintf("Warning message:\n Information matrix using a COX weight is not invertable.\n Program stopped.\n\n");
			conv[0] = 0.;
			return;
		}
	}

	//memory cleanup
	for(int i=0; i<n; i++) delete[] wdn3[i];
	delete[] wdn3;
	for(int i=0; i<n; i++) delete[] wy3[i];
	delete[] wy3;
	for(int i=0; i<n; i++) delete[] dnci[i];
	delete[] dnci;
	for(int i=0; i<n; i++) delete[] s1s0c[i];
	delete[] s1s0c;
	for(int i=0; i<ncovc; i++) delete[] ric[i];
	delete[] ric;	
}


// .C interface for "Stratified-KM"				  
void FGweight_KM_Strata(float *ts, int *ics, float *zs, int *N, int *NCOV, float *gfg1, float *var1, float *tbase_km, float *lambda10km, float *lambda10sd_km, float *Wlambda, float *Wbeta, float *conv, int *NSTRATA, int *strata_var, float *variance, float *var_conservative) {
	//variable declaration
	int n = *N;
	int ncov = *NCOV;
	int nstrata = *NSTRATA;
	int nfgdiv = 0, nit = 0;
	int njp = 0.;
	int *njpp = &njp;	
	float iflag = 0.;
	float *ifl = &iflag;	
	float max_error = 1.;

	vector<float> tjp(n);
	vector<float> g1(ncov);
	vector<float> g2(ncov);
	vector<float> error(ncov);
	
	float **dnc = new float*[n];
	for(int i=0; i<n; i++) dnc[i] = new float[nstrata];	
	float **yc = new float*[n];
	for(int i=0; i<n; i++) yc[i] = new float[nstrata];	
	float **wdn1 = new float*[n];
	for(int i=0; i<n; i++) wdn1[i] = new float[n];
	float **wy1 = new float*[n];
	for(int i=0; i<n; i++) wy1[i] = new float[n];	
	float **dnci = new float*[n];
	for(int i=0; i<n; i++) dnci[i] = new float[n];
	
	FG_KM_Strata W(n,ncov,nstrata,ts,ics,zs,strata_var);
	W.riskst(wdn1, wy1, dnc, yc, dnci, njpp, tjp);
	for(int i=0; i<ncov; i++) {
		g1[i] = 0.;
		g2[i] = 0.;
	}
	nit = 0;
	
	//newton-raphson for FG model
	while(max_error > 0.00001 && nit < 50) {
		W.fgest1(njp, wy1, wdn1, g1, g2, ifl, ncov);
		nit++;
		if(iflag==1) {
			nfgdiv++;
			Rprintf("Algorithm did not converge when fitting the model with a stratified Kaplan-Meier weight\n Program stopped.\n\n");
			conv[0] = 0.;
			return;
		}
		else {
			max_error = fabs(g2[0]-g1[0]);
			for(int k=1; k<ncov; k++) {
				error[k] = fabs(g2[k]-g1[k]);		
				if(error[k] > max_error) max_error = error[k];
			}
			for(int k=0; k<ncov; k++) g1[k] = g2[k];		
		}//endif	
	}//end while
	if (nit==50) {
		nfgdiv++;
		Rprintf("Algorithm did not converge when fitting the model with a stratified Kaplan-Meier weight\n Program stopped.\n\n");
		conv[0] = 0.;
		return;
	}	
	for(int k=0; k<ncov; k++) gfg1[k]=g2[k];	
	
	//variance calculation
	if(variance[0] == 1.) {
		W.fgest2km(tjp,njp,wy1,wdn1,gfg1,dnc,yc,dnci,var1,ncov,lambda10km,tbase_km,lambda10sd_km,Wlambda,Wbeta,var_conservative);
	}
	
	//memory cleanup
	for(int i=0; i<n; i++) delete[] dnc[i];
	delete[] dnc;
	for(int i=0; i<n; i++) delete[] yc[i];
	delete[] yc;
	for(int i=0; i<n; i++) delete[] wdn1[i];
	delete[] wdn1;
	for(int i=0; i<n; i++) delete[] wy1[i];
	delete[] wy1;
	for(int i=0; i<n; i++) delete[] dnci[i];
	delete[] dnci;
}


// .C interface for "Stratified-COX"
void FGweight_COX_Strata(float *ts, int *ics, float *zs, float *zcs, int *N, int *NCOV, int *NCOVC, float *gfg3, float *var3, float *tbase_cox, float *lambda10cox, float *lambda10sd_cox, float *Wlambda, float *Wbeta, float *conv, float *censdet, int *NSTRATA, int *strata_var, float *variance, float *var_conservative) {
	//variable declaration
	int n = *N;
	int ncov = *NCOV;
	int ncovc = *NCOVC;
	int nstrata = *NSTRATA;
	int ncsdiv = 0, nfgdiv = 0, nit = 0;
	int njp = 0.;
	int *njpp = &njp;		

	float iflagc = 0.;
	float *iflc = &iflagc;
	float iflag = 0.;
	float *ifl = &iflag;
	float iflagric = 0.;
	float *iflric = &iflagric;	
	float max_error = 1.;
	
	vector<float> tjp(n);
	vector<float> g1(ncov);
	vector<float> g2(ncov);
	vector<float> betac(ncovc);
	vector<float> error(ncov);

	float **dnc = new float*[n];
	for(int i=0; i<n; i++) dnc[i] = new float[nstrata];	
	float **yc = new float*[n];
	for(int i=0; i<n; i++) yc[i] = new float[nstrata];	
	float **wdn3 = new float*[n];
	for(int i=0; i<n; i++) wdn3[i] = new float[n];	
	float **wy3 = new float*[n];
	for(int i=0; i<n; i++) wy3[i] = new float[n];		
	float **dnci = new float*[n];
	for(int i=0; i<n; i++) dnci[i] = new float[n];
	float **ric = new float*[ncovc];
	for(int i=0; i<ncovc; i++) ric[i] = new float[ncovc];
	float **rlamc0 = new float*[n];
	for(int i=0; i<n; i++) rlamc0[i] = new float[nstrata];

	float ***s1s0c = new float**[n];
	for(int i=0; i<n; i++) s1s0c[i] = new float*[n];
	for(int i=0; i<n; i++) {
		for(int j=0; j<ncovc; j++) s1s0c[i][j] = new float[nstrata];
	}
	
	FG_COX_Strata W(n,ncov,ncovc,nstrata,ts,ics,zs,zcs,strata_var);
	W.riskst(wdn3, wy3, iflc, betac, dnc, yc, dnci, njpp, tjp, ric, rlamc0, s1s0c, ncovc, censdet);
	if(iflagc==1) {
		ncsdiv++;
		Rprintf("Warning message:\n Algorithm did not converge when fitting the censoring distribution.\n Program stopped.\n\n");
		conv[0] = 0.;
		return;
	}
	for(int i=0; i<ncov; i++) {
		g1[i] = 0.;
		g2[i] = 0.;
	}
	nit = 0;

	//newton-raphson for FG model
	while(max_error > 0.00001 && nit < 50) {
		W.fgest1(njp, wy3, wdn3, g1, g2, ifl, ncov);
		nit++;
		if(iflag==1) {
			nfgdiv++;
			Rprintf("Algorithm did not converge when fitting the model with a stratified-Cox weight\n Program stopped.\n\n");
			conv[0] = 0.;
			return;
		}
		else {
			max_error = fabs(g2[0]-g1[0]);
			for(int k=1; k<ncov; k++) {
				error[k] = fabs(g2[k]-g1[k]);		
				if(error[k] > max_error) max_error = error[k];
			}
			for(int k=0; k<ncov; k++) g1[k] = g2[k];		
		}//endif	
	}//end while
	if (nit==50) {
		nfgdiv++;
		Rprintf("Algorithm did not converge when fitting the model with a Cox weight\n Program stopped.\n\n");
		conv[0] = 0.;
		return;
	}	
	for(int k=0; k<ncov; k++) gfg3[k]=g2[k];	

	//variance calculation
	if(variance[0] == 1.) {
		W.fgest2cox(tjp,njp,betac,ric,rlamc0,s1s0c,wy3,wdn3,gfg3,dnc,yc,dnci,var3,lambda10cox,tbase_cox,lambda10sd_cox,Wlambda,Wbeta,iflric,var_conservative);
		if (iflagric==1) {
			Rprintf("Warning message:\n Information matrix using a COX weight is not invertable.\n Program stopped.\n\n");
			conv[0] = 0.;
			return;
		}
	}

	//memory cleanup
	for(int i=0; i<n; i++) delete[] dnc[i];
	delete[] dnc;
	for(int i=0; i<n; i++) delete[] yc[i];
	delete[] yc;
	for(int i=0; i<n; i++) delete[] wdn3[i];
	delete[] wdn3;
	for(int i=0; i<n; i++) delete[] wy3[i];
	delete[] wy3;
	for(int i=0; i<n; i++) delete[] dnci[i];
	delete[] dnci;
	for(int i=0; i<ncovc; i++) delete[] ric[i];
	delete[] ric;	
	for(int i=0; i<n; i++) delete[] rlamc0[i];
	delete[] rlamc0;

	for(int i=0; i<n; i++) {
		for(int j=0; j<ncovc; j++) delete[] s1s0c[i][j];
		delete[] s1s0c[i];
	}
	delete[] s1s0c;
}


