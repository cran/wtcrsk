//p1.h
#ifndef P1_H
#define P1_H
#include "tntjama/tnt.h"
#include <assert.h>
#include "tntjama/tnt_array1d.h"
#include "tntjama/tnt_array2d.h"
#include "tntjama/tnt_array2d_utils.h"
#include "tntjama/jama_lu.h"
#include <vector>

using namespace TNT;
using namespace JAMA;
using namespace std;

template<class T>
TNT::Array2D<T> invert(const TNT::Array2D<T> &M)
{
	assert(M.dim1() == M.dim2()); // square matrices only please
	JAMA::LU<T> lu(M); // solve for inverse with LU decomposition

	// create identity matrix
	TNT::Array2D<T> id(M.dim1(), M.dim2(), (T)0);
	for (int i = 0; i < M.dim1(); i++) id[i][i] = 1;   
	
	return lu.solve(id); // solves A * A_inv = Identity
}

template<class T>
TNT::Array2D<T> transpose(const TNT::Array2D<T> &M)
{
	TNT::Array2D<T> tran(M.dim2(), M.dim1() );
	for(int r=0; r<M.dim1(); ++r)
		for(int c=0; c<M.dim2(); ++c)
			tran[c][r] = M[r][c];
	return tran;
}

class FG_KM{
  int n;       // sample size
  int ncov;       // number of covariates
  
  Array1D<float> ts; //sorted event time (with tied)
  Array1D<int> ics; //sorted cause indices
  Array2D<float> zs; // covariates

public:
  FG_KM(int N, int NCOV, float *TS, int *ICS, float *ZS);
  void riskst(float **wdn1, float **wy1, vector<float> &dnc, vector<float> &yc, float **dnci, int *njp, vector<float> &tjp) const; 
  void fgest1(int njp, float **wy, float **wdn, vector<float> &gfg1, vector<float> &gfg2, float *iflag, int nvar) const;
  void fgest2km(vector<float> &tjp, int njp, float **wy, float **wdn, float gfg1[], vector<float> &dnc, vector<float> &yc, float **dnci, float var[], int nvar, float lambda10[], float tbase[], float lambda10sd[], float Wlambda[], float Wbeta[], float *var_conservative) const;
};

class FG_COX{
  int n;       // sample size
  int ncov;       // number of covariates
  int ncovc;	//number of covariates for censoring distribution
  
  Array1D<float> ts; //sorted event time (with tied)
  Array1D<int> ics; //sorted cause indices
  Array2D<float> zs; // covariates
  Array2D<float> zcs; // covariatees for censoring distribution

public:
  FG_COX(int N, int NCOV, int NCOVC, float *TS, int *ICS, float *ZS, float *ZCS);
  void riskst(float **wdn3, float **wy3, float *iflagc, vector<float> &betac, vector<float> &dnc, vector<float> &yc, float **dnci, int *njp, vector<float> &tjp, float **ric, vector<float> &rlamc0, float **s1s0c, int nvar, float *censdet) const; 
  void fgest1(int njp, float **wy, float **wdn, vector<float> &gfg1, vector<float> &gfg2, float *iflag, int nvar) const;
  void fgest2cox(vector<float> &tjp, int njp, vector<float> &betac, float **ric, vector<float> &rlamc0, float **s1s0c, float **wy, float **wdn, float gfg1[], vector<float> &dnc, vector<float> &yc, float **dnci, float var[], float lambda10[], float tbase[], float lambda10sd[], float Wlambda[], float Wbeta[], float *iflagric, float *var_conservative) const;
};

class FG_KM_Strata{
  int n;       // sample size
  int ncov;       // number of covariates
  int nstrata;		// number of strata
  
  Array1D<float> ts; //sorted event time (with tied)
  Array1D<int> ics; //sorted cause indices
  Array1D<int> strata_var; //sorted strata values
  Array2D<float> zs; // covariates

public:
  FG_KM_Strata(int N, int NCOV, int NSTRATA, float *TS, int *ICS, float *ZS, int *STRATA_VAR);
  void riskst(float **wdn1, float **wy1, float **dnc, float **yc, float **dnci, int *njp, vector<float> &tjp) const; 
  void fgest1(int njp, float **wy, float **wdn, vector<float> &gfg1, vector<float> &gfg2, float *iflag, int nvar) const;
  void fgest2km(vector<float> &tjp, int njp, float **wy, float **wdn, float gfg1[], float **dnc, float **yc, float **dnci, float var[], int nvar, float lambda10[], float tbase[], float lambda10sd[], float Wlambda[], float Wbeta[], float *var_conservative) const;
};

class FG_COX_Strata{
  int n;       // sample size
  int ncov;       // number of covariates
  int ncovc;	//number of covariates for censoring distribution
  int nstrata;		// number of strata
  
  Array1D<float> ts; //sorted event time (with tied)
  Array1D<int> ics; //sorted cause indices
  Array1D<int> strata_var; //sorted strata values
  Array2D<float> zs; // covariates
  Array2D<float> zcs; // covariatees for censoring distribution

public:
  FG_COX_Strata(int N, int NCOV, int NCOVC, int NSTRATA, float *TS, int *ICS, float *ZS, float *ZCS, int *STRATA_VAR);
  void riskst(float **wdn3, float **wy3, float *iflagc, vector<float> &betac, float **dnc, float **yc, float **dnci, int *njp, vector<float> &tjp, float **ric, float **rlamc0, float ***s1s0c, int nvar, float *censdet) const; 
  void fgest1(int njp, float **wy, float **wdn, vector<float> &gfg1, vector<float> &gfg2, float *iflag, int nvar) const;
  void fgest2cox(vector<float> &tjp, int njp, vector<float> &betac, float **ric, float **rlamc0, float ***s1s0c, float **wy, float **wdn, float gfg1[], float **dnc, float **yc, float **dnci, float var[], float lambda10[], float tbase[], float lambda10sd[], float Wlambda[], float Wbeta[], float *iflagric, float *var_conservative) const;
};
#endif




