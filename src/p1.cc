//p1.cc
#include <cstdlib>
#include <cmath>
#include <R.h>
#include <Rmath.h>
#include <new>
#include "p1.h"

using namespace std;
using namespace TNT;
using namespace JAMA;

float Determinant(float **a,int n) {
	int j2;
	float det = 0.;
	float **m = NULL;
	
	if (n<1) {Rprintf("Dimension too small!\n"); } 
	else if (n==1) {det = a[0][0];} 
	else if (n==2) {
		det = a[0][0]*a[1][1]-a[1][0]*a[0][1];
	} else {
		det = 0.;
		for(int j1=0; j1<n; j1++) {
			m = new float*[n-1];
			for(int i=0; i<n-1; i++) m[i] = new float[n-1];
			for(int i=1; i<n; i++) {
				j2 = 0;
				for(int j=0; j<n; j++) {
                if (j==j1) continue;	
				m[i-1][j2] = a[i][j];
                j2++;
            }
         }
         det += pow(-1.0,1.0+j1+1.0)*a[0][j1]*Determinant(m,n-1);
         for(int i=0; i<n-1; i++) free(m[i]);
         free(m);
      }
   }
   return(det);
}

// Code for KM implementation
FG_KM::FG_KM(int N, int NCOV, float *TS, int *ICS, float *ZS) {
	n = N;
	ncov = NCOV;
	ts = Array1D<float> (n,0.);
	ics = Array1D<int> (n,0.);
	zs = Array2D<float> (n,ncov,0.);
	
	for(int i=0; i<n; i++) {
		ts[i] = TS[i];
		ics[i] = ICS[i];
	}
	for(int i=0; i<n; i++) 
		for(int j=0; j<ncov; j++) zs[i][j] = ZS[i*ncov+j];
}

void FG_KM::riskst(float **wdn1, float **wy1, vector<float> &dnc, vector<float> &yc, float **dnci, int *njp, vector<float> &tjp) const{
	//variable and arrays definition
	vector<int> idxt(n);
	vector<float> gkm(n);

	//fit censoring distribution
	tjp[0] = ts[0];
	*njp += 1;
	idxt[0] = *njp;

	for(int i=1; i<n; i++) {
		if (ts[i] > tjp[*njp-1]) {
			tjp[*njp] = ts[i];
			*njp += 1;
		}
		idxt[i] = *njp-1;
	}
	
	//calculate weights based on risk set	
	for(int j=0; j<*njp; j++) {
		dnc[j] = 0.;
		yc[j] = 0.;
		for(int i=0; i<n; i++) {
			dnci[i][j] = 0.;
          	if (ts[i] == tjp[j] && ics[i] == 0) dnci[i][j] = 1.;
			if (ts[i] >= tjp[j]) yc[j] += 1;	
			if (ts[i] == tjp[j] && ics[i] == 0) dnc[j] += 1;
		}
		if (fabs(yc[j]) > 0 && j==0) gkm[j] = 1.-dnc[j]/yc[j];
		else if(abs(yc[j] > 0) && j>0) gkm[j] = gkm[j-1]*(1.-dnc[j]/yc[j]);
		else gkm[j] = gkm[j-1]; //endif
		
		for(int i=0; i<n; i++) {
			wdn1[i][j] = 0.;
			wy1[i][j] = 0.;
			
			if (ts[i] >= tjp[j]) {
				if (gkm[j] > 0) wy1[i][j] = 1.;
				if (ts[i] == tjp[j] && ics[i] == 1) {
					if (gkm[j] > 0) wdn1[i][j] = 1.;
				}//endif
			}//endif
			
			if (ts[i] < tjp[j] && ics[i] > 1) {
				if (gkm[idxt[i]] > 0) wy1[i][j] = gkm[j]/gkm[idxt[i]];
			}//endif	
		}//endi			
	}//end j
	
	return;
}				   			

void FG_KM::fgest1(int njp, float **wy, float **wdn, vector<float> &gfg1, vector<float> &gfg2, float *iflag, int nvar) const {
	//variable declaration
	float s0, btmp, det;

	vector<float> bzi(n);
	vector<float> s1(nvar);
	vector<float> s1s0(nvar);

	Array2D<float> s2(nvar,nvar,0.);
	Array2D<float> s2s0(nvar,nvar,0.);
	Array2D<float> s1s0sq(nvar,nvar,0.);	
	Array2D<float> ubeta(nvar,1,0.);
	Array2D<float> rbeta0(nvar,nvar,0.);
	Array2D<float> rbeta0inv(nvar,nvar,0.);
	Array2D<float> diff(nvar,1,0.);
	
	float **rbeta = new float*[nvar];
	for(int i=0; i<nvar; i++) rbeta[i] = new float[nvar];	
	
	//initialization
	for(int k=0; k<nvar; k++) {
		ubeta[k][0] = 0.;
		s1s0[k] = 0.;
		for(int l=0; l<nvar; l++) {
			rbeta[k][l] = 0.;
			s2s0[k][l] = 0.;
			s1s0sq[k][l] = 0.;
		}
	}//endk

	//calculation
	for (int i=0; i<n; i++) {
		btmp = 0.;	
		for(int k=0; k<nvar; k++) btmp += gfg1[k]*zs[i][k];		
		bzi[i] = exp(btmp);
	}//endi
	
	for(int j=0; j<njp; j++) {
		s0 = 0.;
		for(int k=0; k<nvar; k++) {
			s1[k] = 0.;
			for(int l=0; l<nvar; l++) {
				s2[k][l] = 0.;			
			}
		}//end k

		for (int i=0; i<n; i++) {
			s0 += wy[i][j]*bzi[i];
			for(int k=0; k<nvar; k++) s1[k] += wy[i][j]*zs[i][k]*bzi[i];
			for(int k=0; k<nvar; k++) {
				for(int l=0; l<nvar; l++) s2[k][l] += wy[i][j]*zs[i][k]*zs[i][l]*bzi[i];
			}
		}//end i

		if (fabs(s0) > 0) {
			for(int k=0; k<nvar; k++) s1s0[k] = s1[k]/s0;
          	for(int k=0; k<nvar; k++) {
				for(int l=0; l<nvar; l++) s2s0[k][l] = s2[k][l]/s0;
			}
		}
        else {
			for(int k=0; k<nvar; k++) s1s0[k] = 0.;
          	for(int k=0; k<nvar; k++) {
				for(int l=0; l<nvar; l++) s2s0[k][l] = 0.;
			}
		}//endif

		for(int k=0; k<nvar; k++) {
			for(int l=0; l<nvar; l++) s1s0sq[k][l] = s1s0[k]*s1s0[l];
		}
		
		for(int i=0; i<n; i++) {
			for(int k=0; k<nvar; k++) ubeta[k][0] += (zs[i][k]-s1s0[k])*wdn[i][j]; 
			for(int k=0; k<nvar; k++) {
				for(int l=0; l<nvar; l++) rbeta[k][l] += (s2s0[k][l]-s1s0sq[k][l])*wdn[i][j];
			}
		}//end i
	}//end j	
	
	
	for(int k=0; k<nvar; k++)
		for(int l=0; l<nvar; l++)
			rbeta0[k][l] = rbeta[k][l];

	det = Determinant(rbeta,nvar);
	if(fabs(det) > 0) {
		rbeta0inv = invert(rbeta0);	
		diff =  matmult(rbeta0inv,ubeta);
		for(int k=0; k<nvar; k++) {
			gfg2[k] = gfg1[k]+diff[k][0];
		}
	}
	else {
		*iflag = 1.;
	}			
	
	//memory cleanup
	for(int i=0; i<nvar; i++) delete[] rbeta[i];
	delete[] rbeta;
	
	return;
}

void FG_KM::fgest2km(vector<float> &tjp, int njp, float **wy, float **wdn, float gfg1[], vector<float> &dnc, vector<float> &yc, float **dnci, float var[], int nvar, float lambda10[], float tbase[], float lambda10sd[], float Wlambda[], float Wbeta[]) const {
	//variable declaration	
	int count = 0, count2 = 0, tj = 0, count4 = 0;
	float wdnsum, btmp, event_or_not, det, jnk4, jnk5, jnk6;
	vector<float> bzi(n);
	vector<float> a10hat(njp);
	vector<float> achat(n);
	vector<float> s0(njp);
	vector<float> s1(nvar);
	vector<float> EwdnI(nvar);
	
	Array2D<float> etahat(n,nvar,0.);
	Array2D<float> psihat(n,nvar,0.);
	Array2D<float> s1s0(njp,nvar,0.);
	Array2D<float> temp(n,nvar,0.);
	Array2D<float> tempsum2(n,nvar,0.);
	Array2D<float> Ewdn(njp,nvar,0.);
	
	Array2D<float> s2(nvar,nvar,0.);
	Array2D<float> s2s0(nvar,nvar,0.);
	Array2D<float> s1s0sq(nvar,nvar,0.);
	Array2D<float> vcov(nvar,nvar,0.);
	Array2D<float> tempsum(nvar,nvar,0.);
	
	Array2D<float> omega0(nvar,nvar,0.);
	Array2D<float> omega0inv(nvar,nvar,0.);
	Array2D<float> sigma(nvar,nvar,0.);
	Array2D<float> final(nvar,nvar,0.);
	
	Array3D<float> temp2(n,nvar,nvar,0.);
	Array3D<float> omegahat(n,nvar,nvar,0.);
	
	float **wdm1 = new float*[n];
	for(int i=0; i<n; i++) wdm1[i] = new float[njp];	
	float **yci = new float*[n];
	for(int i=0; i<n; i++) yci[i] = new float[njp];	
	float **dmci = new float*[n];
	for(int i=0; i<n; i++) dmci[i] = new float[njp];	
	float **W_lambda = new float*[n];
	for(int i=0; i<n; i++) W_lambda[i] = new float[njp];	
	float **W_lambda1 = new float*[n];
	for(int i=0; i<n; i++) W_lambda1[i] = new float[njp];
	float **W_lambda2 = new float*[n];
	for(int i=0; i<n; i++) W_lambda2[i] = new float[njp];
	float **Wc = new float*[n];
	for(int i=0; i<n; i++) Wc[i] = new float[njp];
	float **omegasum = new float*[nvar];
	for(int i=0; i<nvar; i++) omegasum[i] = new float[nvar];
	
	//initializaztion
	for(int i=0; i<n; i++) {
		btmp = 0.;	
		for(int k=0; k<nvar; k++) btmp += gfg1[k]*zs[i][k];		
		bzi[i] = exp(btmp);
		for(int k=0; k<nvar; k++) {
			etahat[i][k] = 0.;
			psihat[i][k] = 0.;
			for(int l=0; l<nvar; l++) omegahat[i][k][l]=0.;
		}
	}	
	
	for(int i=0; i<n; i++) {
		for(int k=0; k<nvar; k++) {
			for(int l=0; l<nvar; l++) temp2[i][k][l] = 0.;	
		}
	}
	
	for(int k=0; k<nvar; k++) 
		for(int l=0; l<nvar; l++) omegasum[k][l] = 0.;	

	for(int i=0; i<n; i++) {
		for(int k=0; k<nvar; k++) tempsum2[i][k] = 0.;
	}		
	
	for(int i=0; i<n; i++) {
		for(int k=0; k<nvar; k++) {
			for(int l=0; l<nvar; l++) tempsum[k][l] = 0.;
		}
	}
	
	for(int k=0; k<nvar; k++) EwdnI[k] = 0.;
	
	//calculation
	for(int j=0; j<njp; j++) {
		s0[j] = 0.;
		for(int k=0; k<nvar; k++) {
			s1[k] = 0.;
			for(int l=0; l<nvar; l++) {
				s2[k][l] = 0.;			
			}
		}//end k,l
		wdnsum = 0.;	
		for(int i=0; i<n; i++) {
			s0[j] += wy[i][j]*bzi[i];
			for(int k=0; k<nvar; k++) s1[k] += wy[i][j]*zs[i][k]*bzi[i];
			for(int k=0; k<nvar; k++) {
				for(int l=0; l<nvar; l++) s2[k][l] += wy[i][j]*zs[i][k]*zs[i][l]*bzi[i];
			}
			wdnsum += wdn[i][j];
		}//endi
		
		if(fabs(s0[j]) > 0) {	
			for(int k=0; k<nvar; k++) s1s0[j][k] = s1[k]/s0[j];
          	for(int k=0; k<nvar; k++) {
				for(int l=0; l<nvar; l++) s2s0[k][l] = s2[k][l]/s0[j];
			}
			if(j > 0) a10hat[j] = a10hat[j-1]+wdnsum/s0[j];
			else a10hat[j] = wdnsum/s0[j];
		}
		else {
			for(int k=0; k<nvar; k++) s1s0[j][k] = 0.;
          	for(int k=0; k<nvar; k++) {
				for(int l=0; l<nvar; l++) s2s0[k][l] = 0.;
			}	
			if(j > 0) a10hat[j] = a10hat[j-1];
			else a10hat[j] = 0.;
		}//endif
	
		//estimation of baseline hazard
		event_or_not = 0.;
		for(int i=0; i<n; i++) {
			if(ics[i]==1 && wdn[i][j]>0) {
				event_or_not = 1.;
				break;
			}
		}
		if(event_or_not == 1.) {
			tbase[count] = tjp[j];
			if(count > 0 && s0[j] > 0) lambda10[count] = lambda10[count-1]+wdnsum/s0[j];
			else if(s0[j] > 0) lambda10[count] = wdnsum/s0[j];
			count++;
		}//endif
		
		for(int k=0; k<nvar; k++) {
			for(int l=0; l<nvar; l++) s1s0sq[k][l] = s1s0[j][k]*s1s0[j][l];
		}
		
		for(int i=0; i<n; i++) {
			if(j > 0) wdm1[i][j] = wdn[i][j]-wy[i][j]*bzi[i]*(a10hat[j]-a10hat[j-1]);
			else wdm1[i][j] = wdn[i][j]-wy[i][j]*bzi[i]*a10hat[j];
			for(int k=0; k<nvar; k++) etahat[i][k] += (zs[i][k]-s1s0[j][k])*wdm1[i][j];
			for(int k=0; k<nvar; k++) {
				for(int l=0; l<nvar; l++) omegahat[i][k][l] += (s2s0[k][l]-s1s0sq[k][l])*wdn[i][j];
			}//end k,l
		}//endi			
		
		if(event_or_not == 1.) {         //calculate W_lambda1
			for(int i=0; i<n; i++) {
				if(count2 > 0 && s0[j] > 0) W_lambda1[i][count2] = W_lambda1[i][count2-1]+wdm1[i][j]/s0[j];
				else if( count2 == 0 && s0[j] > 0) W_lambda1[i][count2] = wdm1[i][j]/s0[j];
			}
			count2++;
		}	
	}//endj
	
	//calculate dmci & yci
	for(int j=0; j<njp; j++) {
		if(j > 0) achat[j] = achat[j-1]+dnc[j]/yc[j];
		else achat[j] = dnc[j]/yc[j];
		for(int i=0; i<n; i++) {
			yci[i][j] = 0.;
			if (ts[i] >= tjp[j]) yci[i][j] = 1.;
			dmci[i][j] = dnci[i][j]-yci[i][j]*(achat[j]-achat[j-1]);
		}
	}//endj
	
	//calculate Wc matrix and variance
	for(int j=0; j<njp; j++) {
		for(int i=0; i<n; i++) {
			if(j > 0) Wc[i][j] = Wc[i][j-1]+dmci[i][j]/yc[j];
			else Wc[i][j] = dmci[i][j]/yc[j];
		}
	}//endj	
	for(int i=0; i<n; i++) {
		for(int j=0; j<n; j++) {
			for(int k=0; k<nvar; k++) temp[j][k] = 0.;
			for(int j2=0; j2<njp; j2++) {
				if(ts[j] == tjp[j2]) {
					tj = j2;
					break;
				}
			}//end j2
			for(int j2=0; j2<njp; j2++) {
				if(tj < j2) {
					for(int k=0; k<nvar; k++) temp[j][k] += (zs[j][k]-s1s0[j2][k])*wdm1[j][j2]*(Wc[i][tj]-Wc[i][j2]);
				}
			}//end j2
			for(int k=0; k<nvar; k++) psihat[i][k] += temp[j][k];
		}
	}//end i,j

	for(int i=0; i<n; i++) {
		for(int k=0; k<nvar; k++) {
			tempsum2[i][k] = etahat[i][k]+psihat[i][k];                      
 			for(int l=0; l<nvar; l++) {
				temp2[i][k][l] += (etahat[i][k]+psihat[i][k])*(etahat[i][l]+psihat[i][l]);
				omegasum[k][l] += omegahat[i][k][l];
			}
		} 	
	}//endi	
	
	for(int i=0; i<n; i++) {
		for(int k=0; k<nvar; k++) {
			for(int l=0; l<nvar; l++) tempsum[k][l] += temp2[i][k][l];
		}
	}	

	det = Determinant(omegasum,nvar);
	for(int k=0; k<nvar; k++) {
		for(int l=0; l<nvar; l++) {
			omega0[k][l] = omegasum[k][l];
			sigma[k][l] = tempsum[k][l];
		}
	}
	if(fabs(det) > 0) {
		omega0inv = invert(omega0);	
		final = matmult(matmult(omega0inv,sigma),omega0inv);
		for(int k=0; k<nvar; k++) {
			for(int l=0; l<nvar; l++) vcov[k][l] = final[k][l];
		}
		for(int k=0; k<nvar; k++) var[k] = vcov[k][k];		
	}
	else {
		Rprintf("Singular matrix!\n");
		return;
	}
	
	//variance for baseline hazard
	for(int j=0; j<njp; j++) {
		event_or_not = 0.;
		for(int i=0; i<n; i++) {
			if(ics[i]==1 && wdn[i][j]>0) {
				event_or_not = 1.;
				break;
			}
		}
		if(event_or_not ==1.) {
			for(int k=0; k<nvar; k++) {
				if(count4 > 0) Ewdn[count4][k] = Ewdn[count4-1][k];
				else if(count4 == 0) Ewdn[count4][k] = 0.;
			}
			for(int i=0; i<n; i++) {
				if(s0[j] > 0) {
					for(int k=0; k<nvar; k++) Ewdn[count4][k] += s1s0[j][k]*wdn[i][j]/s0[j];
				}
			}
			
			for(int k=0; k<nvar; k++) EwdnI[k] = 0.;
			for(int k=0; k<nvar; k++) {
				for(int l=0; l<nvar; l++) EwdnI[k] += Ewdn[count4][l]*omega0inv[l][k];
			}

			for(int i=0; i<n; i++) {	
				jnk4 = 0.;
				for(int k=0; k<nvar; k++) jnk4 += EwdnI[k]*tempsum2[i][k];
				W_lambda2[i][count4] = jnk4;
			}	
			count4++;			
		}
	}//endj	
	
	for(int i=0; i<n; i++) {
		for(int j=0; j<count4; j++) W_lambda[i][j] = W_lambda1[i][j]-W_lambda2[i][j];
	}
	for(int i=0; i<n; i++) {
		for(int j=0; j<count4; j++) Wlambda[i*count4+j] = W_lambda[i][j];
	}	
	
	for(int i=0; i<n; i++) {
		for(int k=0; k<nvar; k++) {
			jnk6 = 0.;
			for (int l=0; l<nvar; l++) jnk6 += omega0inv[k][l]*tempsum2[i][l];
			Wbeta[i*nvar+k] = jnk6;
		}		
	}
	for(int j=0; j<count4; j++) {
		jnk5 = 0.;
		for(int i=0; i<n; i++) jnk5 += pow((double)W_lambda[i][j],2.);
		lambda10sd[j] = sqrt(jnk5);
	}	
	
	//memory cleanup
	for(int i=0; i<n; i++) delete[] wdm1[i];
	delete[] wdm1;
	for(int i=0; i<n; i++) delete[] yci[i];
	delete[] yci;
	for(int i=0; i<n; i++) delete[] dmci[i];
	delete[] dmci;
	for(int i=0; i<n; i++) delete[] W_lambda[i];
	delete[] W_lambda;
	for(int i=0; i<n; i++) delete[] W_lambda1[i];
	delete[] W_lambda1;
	for(int i=0; i<n; i++) delete[] W_lambda2[i];
	delete[] W_lambda2;
	for(int i=0; i<n; i++) delete[] Wc[i];
	delete[] Wc;
	for(int i=0; i<nvar; i++) delete[] omegasum[i];
	delete[] omegasum;
}


// Code for Cox implementation
FG_COX::FG_COX(int N, int NCOV, int NCOVC, float *TS, int *ICS, float *ZS, float *ZCS) {
	n = N;
	ncov = NCOV;
	ncovc = NCOVC;
	ts = Array1D<float> (n,0.);
	ics = Array1D<int> (n,0.);
	zs = Array2D<float> (n,ncov,0.);
	zcs = Array2D<float> (n,ncovc,0.);
	
	for(int i=0; i<n; i++) {
		ts[i] = TS[i];
		ics[i] = ICS[i];
	}
	for(int i=0; i<n; i++) {
		for(int j=0; j<ncov; j++) zs[i][j] = ZS[i*ncov+j];
		for(int j=0; j<ncovc; j++) zcs[i][j] = ZCS[i*ncovc+j];
	}
}

void FG_COX::riskst(float **wdn3, float **wy3, float *iflagc, vector<float> &betac, vector<float> &dnc, vector<float> &yc, float **dnci, int *njp, vector<float> &tjp, float **ric, vector<float> &rlamc0, float **s1s0c, int nvar, float *censdet) const{
	//variable and arrays definition
	int nit = 0;
	float s0, tmp, max_error=1.;
	
	vector<int> idxt(n);
	vector<float> bzi(n);
	vector<float> error(nvar);
	vector<float> s1(nvar);
	vector<float> bc1(nvar);
	vector<float> bc2(nvar);
	
	Array2D<float> s2(nvar,nvar,0.);
	Array2D<float> s1s0sq(nvar,nvar,0.);
	Array2D<float> s2s0(nvar,nvar,0.);
	
	void ccox(int n, const Array1D<float> &ts, const Array1D<int> &ics, const Array2D<float> &zs, vector<float> &tjp, int njp, vector<float> &bc1, vector<float> &bc2, float **dnci, float *iflagc, int nvar, float *censdet);
	
	//dynamic arrays allocation
	float **gsurv = new float*[n];
	for(int i=0; i<n; i++) gsurv[i] = new float[n];

	//fit censoring distribution
	for(int k=0; k<nvar; k++) {
		bc1[k] = 0.;
		bc2[k] = 0.;
	}

	tjp[0] = ts[0];
	*njp += 1;
	idxt[0] = *njp;

	for(int i=1; i<n; i++) {
		if (ts[i] > tjp[*njp-1]) {
			tjp[*njp] = ts[i];
			*njp += 1;
		}
		idxt[i] = *njp-1;
	}
	
	while(max_error > 0.00001 && nit < 50) {
		ccox(n, ts, ics, zcs, tjp, *njp, bc1, bc2, dnci, iflagc, nvar, censdet);
		nit++;
		if(*iflagc==1) return;
		else {
			max_error = fabs(bc2[0]-bc1[0]);
			for(int k=1; k<nvar; k++) {
				error[k] = fabs(bc2[k]-bc1[k]);		
				if(error[k] > max_error) max_error = error[k];
			}
			for(int k=0; k<nvar; k++) bc1[k] = bc2[k];		
		}//endif	
	}//end while
	
	if (nit==50) {
		*iflagc = 1;
		return;
	}	
	for(int k=0; k<nvar; k++) betac[k]=bc2[k];	
	
	//calculate weights based on risk set
	for (int i=0; i<n; i++) {
		tmp = 0.;	
		for(int k=0; k<nvar; k++) {
			tmp += betac[k]*zcs[i][k];
		}
		bzi[i] = exp(tmp);
	}
	for(int k=0; k<nvar; k++)
		for(int l=0; l<nvar; l++) ric[k][l] = 0.;	
	
	for(int j=0; j<*njp; j++) {
		dnc[j] = 0.;
		yc[j] = 0.;
		s0 = 0.;
		for(int k=0; k<nvar; k++) {
			s1[k] = 0.;
			for(int l=0; l<nvar; l++) {
				s2[k][l] = 0.;			
			}
		}
		for(int i=0; i<n; i++) {
			if (ts[i] >= tjp[j]) {
				yc[j] += 1;
				s0 += bzi[i];
				for(int k=0; k<nvar; k++) s1[k] += zcs[i][k]*bzi[i];
				for(int k=0; k<nvar; k++) {
					for(int l=0; l<nvar; l++) s2[k][l] += zcs[i][k]*zcs[i][l]*bzi[i]; //****
				}
			}//endif	
			if (ts[i] == tjp[j] && ics[i] == 0) dnc[j] += 1;
		}
		
		if (fabs(s0) > 0) {
			for(int k=0; k<nvar; k++) s1s0c[j][k] = s1[k]/s0;
          	for(int k=0; k<nvar; k++) {
				for(int l=0; l<nvar; l++) s2s0[k][l] = s2[k][l]/s0;
			}
			if(j == 0) rlamc0[j] = dnc[j]/s0;
			else rlamc0[j] = rlamc0[j-1] + dnc[j]/s0;
		}
        else {
			for(int k=0; k<nvar; k++) s1s0c[j][k] = 0.;
          	for(int k=0; k<nvar; k++) {
				for(int l=0; l<nvar; l++) s2s0[k][l] = 0.;
			}
			if(j > 0) rlamc0[j] = rlamc0[j-1];
		} //endif

		for(int k=0; k<nvar; k++) {
			for(int l=0; l<nvar; l++) s1s0sq[k][l] = s1s0c[j][k]*s1s0c[j][l];
		}		
		for(int k=0; k<nvar; k++) {
			for(int l=0; l<nvar; l++) ric[k][l] += (s2s0[k][l] - s1s0sq[k][l])*dnc[j];
		}
		
		for(int i=0; i<n; i++) {
			gsurv[i][j] = exp(0-rlamc0[j]*bzi[i]);
			wdn3[i][j] = 0.;
			wy3[i][j] = 0.;
			if (ts[i] >= tjp[j]) {
				if (gsurv[i][j] > 0) wy3[i][j] = 1.;
				if (ts[i] == tjp[j] && ics[i] == 1) {
					if (gsurv[i][j] > 0) wdn3[i][j] = 1.;
				}//endif
			}//endif
			if (ts[i] < tjp[j] && ics[i] > 1) {
				if (gsurv[i][idxt[i]] > 0) wy3[i][j] = gsurv[i][j]/gsurv[i][idxt[i]];
			}//endif	
		}//endi			
	}//end j
	
	//memory cleanup
	for(int i=0; i<n; i++) delete[] gsurv[i];
	delete[] gsurv;	
	
	return;
}		

void ccox(int n, const Array1D<float> &ts, const Array1D<int> &ics, const Array2D<float> &zs, vector<float> &tjp, int njp, vector<float> &bc1, vector<float> &bc2, float **dnci, float *iflagc, int nvar, float *censdet) {
	//variable declaration
	float s0, btmp, det;

	vector<float> s1(nvar);
	vector<float> s1s0(nvar);
	
	Array2D<float> s2(nvar,nvar,0.);
	Array2D<float> s2s0(nvar,nvar,0.);
	Array2D<float> s1s0sq(nvar,nvar,0.);
	
	Array2D<float> uc(nvar,1,0.);
	Array2D<float> rc0(nvar,nvar,0.);
	Array2D<float> rc0inv(nvar,nvar,0.);
	Array2D<float> diff(nvar,1,0.);
	
	float **rc = new float*[nvar];
	for(int i=0; i<nvar; i++) rc[i] = new float[nvar];
	
	//initialization
	for(int k=0; k<nvar; k++) {
		uc[k][0] = 0.;
		s1s0[k] = 0.;
		for(int l=0; l<nvar; l++) {
			rc[k][l] = 0.;
			s2s0[k][l] = 0.;
			s1s0sq[k][l] = 0.;
		}
	}
	
	//calculation
	for (int j=0; j<njp; j++) {
		s0 = 0.;
		for(int k=0; k<nvar; k++) {
			s1[k] = 0.;
			for(int l=0; l<nvar; l++) {
				s2[k][l] = 0.;			
			}
		}//end k
		for (int i=0; i<n; i++) {
			btmp = 0.;
			double tmp = 0.;	
			for(int k=0; k<nvar; k++) tmp += bc1[k]*zs[i][k];
			btmp = exp(tmp);

			if (ts[i] >= tjp[j]) {
            	s0 += btmp;
				for(int k=0; k<nvar; k++) s1[k] += zs[i][k]*btmp;
            	for(int k=0; k<nvar; k++) {
					for(int l=0; l<nvar; l++) s2[k][l] += zs[i][k]*zs[i][l]*btmp;
				}
			}//endif
		}//endi
		
	    if (fabs(s0) > 0.) {
			for(int k=0; k<nvar; k++) s1s0[k] = s1[k]/s0;
          	for(int k=0; k<nvar; k++)
				for(int l=0; l<nvar; l++)
          				s2s0[k][l] = s2[k][l]/s0;
		}
        else {
			for(int k=0; k<nvar; k++) s1s0[k] = 0.;
          	for(int k=0; k<nvar; k++)
				for(int l=0; l<nvar; l++)
          				s2s0[k][l] = 0.;
		} //endif
		
 		for(int k=0; k<nvar; k++) {
			for(int l=0; l<nvar; l++) s1s0sq[k][l] = s1s0[k]*s1s0[l];
		}
    
		for (int i=0; i<n; i++) {
			dnci[i][j] = 0.;
          	if (ts[i] == tjp[j] && ics[i] == 0) dnci[i][j] = 1.;
			for(int k=0; k<nvar; k++) {
          			uc[k][0] += (zs[i][k] - s1s0[k])*dnci[i][j];	
			}
			for(int k=0; k<nvar; k++) {
				for(int l=0; l<nvar; l++) {
          				rc[k][l] += (s2s0[k][l] - s1s0sq[k][l])*dnci[i][j];
				}
			}
		}//endi	
	}//end j
	
	for(int k=0; k<nvar; k++)
		for(int l=0; l<nvar; l++)
			rc0[k][l] = rc[k][l];	

	det = Determinant(rc,nvar);
	if(fabs(det) > 0) {
		if(fabs(det) < 0.00001) {
			censdet[0] = 0.;
			iflagc[0] = 1.;
			return;
		}
		rc0inv = invert(rc0);	
		diff =  matmult(rc0inv,uc);
		for(int k=0; k<nvar; k++) {
			bc2[k] = bc1[k]+diff[k][0];
		}
	}
	else {
		*iflagc = 1.;
		return;
	}
	
	//memory cleanup
	for(int i=0; i<nvar; i++) delete[] rc[i];
	delete[] rc;
}

void FG_COX::fgest1(int njp, float **wy, float **wdn, vector<float> &gfg1, vector<float> &gfg2, float *iflag, int nvar) const {
	//variable declaration
	float s0, btmp, det;
	
	vector<float> bzi(n);
	vector<float> s1(nvar);
	vector<float> s1s0(nvar);
	
	Array2D<float> s2(nvar,nvar,0.);
	Array2D<float> s2s0(nvar,nvar,0.);
	Array2D<float> s1s0sq(nvar,nvar,0.);
	
	Array2D<float> ubeta(nvar,1,0.);
	Array2D<float> rbeta0(nvar,nvar,0.);
	Array2D<float> rbeta0inv(nvar,nvar,0.);
	Array2D<float> diff(nvar,1,0.);
	
	float **rbeta = new float*[nvar];
	for(int i=0; i<nvar; i++) rbeta[i] = new float[nvar];	
	
	//initialization
	for(int k=0; k<nvar; k++) {
		ubeta[k][0] = 0.;
		s1s0[k] = 0.;
		for(int l=0; l<nvar; l++) {
			rbeta[k][l] = 0.;
			s2s0[k][l] = 0.;
			s1s0sq[k][l] = 0.;
		}
	}//endk

	//calculation
	for (int i=0; i<n; i++) {
		btmp = 0.;	
		for(int k=0; k<nvar; k++) btmp += gfg1[k]*zs[i][k];		
		bzi[i] = exp(btmp);
	}//endi
	
	for(int j=0; j<njp; j++) {
		s0 = 0.;
		for(int k=0; k<nvar; k++) {
			s1[k] = 0.;
			for(int l=0; l<nvar; l++) {
				s2[k][l] = 0.;			
			}
		}//end k

		for (int i=0; i<n; i++) {
			s0 += wy[i][j]*bzi[i];
			for(int k=0; k<nvar; k++) s1[k] += wy[i][j]*zs[i][k]*bzi[i];
			for(int k=0; k<nvar; k++) {
				for(int l=0; l<nvar; l++) s2[k][l] += wy[i][j]*zs[i][k]*zs[i][l]*bzi[i];
			}
		}//end i

		if (fabs(s0) > 0) {
			for(int k=0; k<nvar; k++) s1s0[k] = s1[k]/s0;
          	for(int k=0; k<nvar; k++) {
				for(int l=0; l<nvar; l++) s2s0[k][l] = s2[k][l]/s0;
			}
		}
        else {
			for(int k=0; k<nvar; k++) s1s0[k] = 0.;
          	for(int k=0; k<nvar; k++) {
				for(int l=0; l<nvar; l++) s2s0[k][l] = 0.;
			}
		}//endif

		for(int k=0; k<nvar; k++) {
			for(int l=0; l<nvar; l++) s1s0sq[k][l] = s1s0[k]*s1s0[l];
		}
		
		for(int i=0; i<n; i++) {
			for(int k=0; k<nvar; k++) ubeta[k][0] += (zs[i][k]-s1s0[k])*wdn[i][j]; 
			for(int k=0; k<nvar; k++) {
				for(int l=0; l<nvar; l++) rbeta[k][l] += (s2s0[k][l]-s1s0sq[k][l])*wdn[i][j];
			}
		}//end i
	}//end j	

	for(int k=0; k<nvar; k++)
		for(int l=0; l<nvar; l++)
			rbeta0[k][l] = rbeta[k][l];

	det = Determinant(rbeta,nvar);
	if(fabs(det) > 0) {
		rbeta0inv = invert(rbeta0);	
		diff =  matmult(rbeta0inv,ubeta);
		for(int k=0; k<nvar; k++) {
			gfg2[k] = gfg1[k]+diff[k][0];
		}
	}
	else {
		*iflag = 1.;
	}			
	
	//memory cleanup
	for(int i=0; i<nvar; i++) delete[] rbeta[i];
	delete[] rbeta;
	
	return;
}

void FG_COX::fgest2cox(vector<float> &tjp, int njp, vector<float> &betac, float **ric, vector<float> &rlamc0, float **s1s0c, float **wy, float **wdn, float gfg1[], vector<float> &dnc, vector<float> &yc, float **dnci, float var[], float lambda10[], float tbase[], float lambda10sd[], float Wlambda[], float Wbeta[], float *iflagric) const {
	//variable declaration
	int count = 0, count2 = 0, count4 = 0, tj;
	float wdnsum, btmp, btmp2, event_or_not, det, Wc1 = 0., Wc2 = 0., jnk1 = 0., jnk4, jnk5, jnk6;
	
	vector<float> bzi(n);
	vector<float> bzic(n);
	vector<float> a10hat(n);
	vector<float> s0(njp);
	vector<float> s1(ncov);
	vector<float> s0c(njp);
	vector<float> EwdnI(ncov);
	
	Array2D<float> s1s0(njp,ncov,0.);
	Array2D<float> s2(ncov,ncov,0.);
	Array2D<float> s2s0(ncov,ncov,0.);
	Array2D<float> s1s0sq(ncov,ncov,0.);
	Array2D<float> temp(njp,ncov,0.);
	Array2D<float> Ewdn(njp,ncov,0.);
	Array2D<float> etahat(n,ncov,0.);
	Array2D<float> psihat(n,ncov,0.);
	Array2D<float> tmp_1(n,ncovc,0.);
	Array2D<float> tempsum(ncov,ncov,0.);
	Array2D<float> tempsum2(n,ncov,0.);
	Array2D<float> vcov(ncov,ncov,0.);
	
	Array2D<float> ric0(ncovc,ncovc,0.);
	Array2D<float> ric0inv(ncovc,ncovc,0.);
	Array2D<float> ht(ncovc,1,0.);
	Array2D<float> jnktmp(1,ncovc,0.);
	Array2D<float> omega0(ncov,ncov,0.);
	Array2D<float> omega0inv(ncov,ncov,0.);
	Array2D<float> sigma(ncov,ncov,0.);
	Array2D<float> final(ncov,ncov,0.);
	
	Array3D<float> temp2(n,ncov,ncov,0.);
	Array3D<float> omegahat(n,ncov,ncov,0.);	
	
	float **Wc = new float*[n];
	for(int i=0; i<n; i++) Wc[i] = new float[njp];	
	float **tmp_2 = new float*[n];
	for(int i=0; i<n; i++) tmp_2[i] = new float[njp];
	float **W_lambda = new float*[n];
	for(int i=0; i<n; i++) W_lambda[i] = new float[njp];
	float **W_lambda1 = new float*[n];
	for(int i=0; i<n; i++) W_lambda1[i] = new float[njp];
	float **W_lambda2 = new float*[n];
	for(int i=0; i<n; i++) W_lambda2[i] = new float[njp];
	float **wdm1 = new float*[n];
	for(int i=0; i<n; i++) wdm1[i] = new float[njp];	
	float **yci = new float*[n];
	for(int i=0; i<n; i++) yci[i] = new float[njp];	
	float **dmci = new float*[n];
	for(int i=0; i<n; i++) dmci[i] = new float[njp];	
	float **omegasum = new float*[ncov];
	for(int i=0; i<ncov; i++) omegasum[i] = new float[ncov];
	
	float ***h = new float**[n];
	for(int i=0; i<n; i++) h[i] = new float*[njp];
	for(int i=0; i<n; i++) {
		for(int j=0; j<njp; j++) h[i][j] = new float[ncovc];
	}
	float ***jnk0 = new float**[n];
	for(int i=0; i<n; i++) jnk0[i] = new float*[njp];
	for(int i=0; i<n; i++) {
		for(int j=0; j<njp; j++) jnk0[i][j] = new float[ncovc];
	}

	//basic calculation and initialization
	for(int i=0; i<n; i++) {
		btmp = 0.;
		btmp2 = 0.;	
		for(int k=0; k<ncov; k++) btmp += gfg1[k]*zs[i][k];
		for(int k=0; k<ncovc; k++) btmp2 += betac[k]*zcs[i][k];

		bzi[i] = exp(btmp);
		bzic[i] = exp(btmp2);
		for(int k=0; k<ncov; k++) {
			etahat[i][k] = 0.;
			psihat[i][k] = 0.;
			for(int l=0; l<ncov; l++) omegahat[i][k][l]=0.;
		}
	}//endi	
	
	for(int i=0; i<n; i++) {
		for(int k=0; k<ncovc; k++) tmp_1[i][k] = 0.;
	}	

	for(int i=0; i<n; i++) {
		for(int k=0; k<ncov; k++) {
			for(int l=0; l<ncov; l++) temp2[i][k][l] = 0.;
		}
	}

	for(int k=0; k<ncov; k++) 
		for(int l=0; l<ncov; l++) omegasum[k][l] = 0.;

	for(int i=0; i<n; i++)
		for(int k=0; k<ncov; k++) tempsum2[i][k] = 0.;

	for(int k=0; k<ncov; k++)
		for(int l=0; l<ncov; l++) tempsum[k][l] = 0.;
	
	for(int k=0; k<ncov; k++) EwdnI[k] = 0.;	
	
	//calculation for variance
	for(int j=0; j<njp; j++) {
		s0[j] = 0.;
		for(int k=0; k<ncov; k++) {
			s1[k] = 0.;
			for(int l=0; l<ncov; l++) s2[k][l] = 0.;			
		}
		wdnsum = 0.;

		for(int i=0; i<n; i++) {
			s0[j] += wy[i][j]*bzi[i];
			for(int k=0; k<ncov; k++) s1[k] += wy[i][j]*zs[i][k]*bzi[i];
			for(int k=0; k<ncov; k++) {
				for(int l=0; l<ncov; l++) s2[k][l] += wy[i][j]*zs[i][k]*zs[i][l]*bzi[i];
			}
			wdnsum += wdn[i][j];
		}//endi	
		
		if (fabs(s0[j]) > 0) {	
			for(int k=0; k<ncov; k++) s1s0[j][k] = s1[k]/s0[j];
          	for(int k=0; k<ncov; k++) {
				for(int l=0; l<ncov; l++) s2s0[k][l] = s2[k][l]/s0[j];
			}
			if(j > 0) a10hat[j] = a10hat[j-1]+wdnsum/s0[j];
			else a10hat[j] = wdnsum/s0[j];
		}
		else {
			for(int k=0; k<ncov; k++) s1s0[j][k] = 0.;
          	for(int k=0; k<ncov; k++) {
				for(int l=0; l<ncov; l++) s2s0[k][l] = 0.;
			}	
			if(j > 0) a10hat[j] = a10hat[j-1];
			else a10hat[j] = 0.;
		}//endif

		event_or_not = 0.;
		for(int i=0; i<n; i++) {
			if(ics[i]==1 && wdn[i][j]>0) {
				event_or_not = 1.;
				break;
			}
		}		
		if(event_or_not == 1) {
			tbase[count] = tjp[j];
			if(count > 0) lambda10[count] = lambda10[count-1]+wdnsum/s0[j];
			else lambda10[count] = wdnsum/s0[j];
			count++;
		}//endif

		for(int k=0; k<ncov; k++) {
			for(int l=0; l<ncov; l++) s1s0sq[k][l] = s1s0[j][k]*s1s0[j][l];
		}

		for(int i=0; i<n; i++) {
			if(j > 0) wdm1[i][j] = wdn[i][j]-wy[i][j]*bzi[i]*(a10hat[j]-a10hat[j-1]);
			else wdm1[i][j] = wdn[i][j]-wy[i][j]*bzi[i]*a10hat[j];
			for(int k=0; k<ncov; k++) etahat[i][k] += (zs[i][k]-s1s0[j][k])*wdm1[i][j];
			for(int k=0; k<ncov; k++) {
				for(int l=0; l<ncov; l++) omegahat[i][k][l] += (s2s0[k][l]-s1s0sq[k][l])*wdn[i][j];
			}
		}//endi	
		
		if(event_or_not==1.) {
			for(int i=0; i<n; i++) {
				if(count2 > 0) W_lambda1[i][count2] = W_lambda1[i][count2-1]+wdm1[i][j]/s0[j];
				else W_lambda1[i][count2] = wdm1[i][j]/s0[j];
			}
			count2++;
		}				
	}//endj

	for(int j=0; j<njp; j++) {
		for(int i=0; i<n; i++) {
			yci[i][j] = 0.;
			if (ts[i] >= tjp[j]) yci[i][j] = 1.;
			if(j > 0) dmci[i][j] = dnci[i][j]-yci[i][j]*bzic[i]*(rlamc0[j]-rlamc0[j-1]);
			else dmci[i][j] = dnci[i][j]-yci[i][j]*bzic[i]*rlamc0[j];
		}
	}//endj

	for(int j=0; j<njp; j++) {
		s0c[j] = 0.;	
		for(int i=0; i<n; i++) {
			for(int k=0; k<ncovc; k++) {
				if(j > 0) h[i][j][k] = h[i][j-1][k]+bzic[i]*(zcs[i][k]-s1s0c[j][k])*(rlamc0[j]-rlamc0[j-1]);
				else h[i][j][k] = bzic[i]*(zcs[i][k]-s1s0c[j][k])*rlamc0[j];
			}
			if(ts[i] >= tjp[j]) s0c[j] += bzic[i];
		}
	}//end j
	
	det = Determinant(ric,ncovc);
	for(int k=0; k<ncovc; k++) {
		for(int l=0; l<ncovc; l++) ric0[k][l] = ric[k][l];
	}
	
	if(fabs(det) > 0) ric0inv = invert(ric0);	
	else {
		iflagric[0] = 1.;
		return;
	}

	for(int i=0; i<n; i++) {		
		for(int j=0; j<njp; j++) {
			for(int k=0; k<ncovc; k++) tmp_1[i][k] += (zcs[i][k]-s1s0c[j][k])*dmci[i][j];
		}			
	}//end i
	
	for(int j=0; j<njp; j++) {
		if(s0c[j] > 0.) {
			for(int i=0; i<n; i++) {
				if(j > 0) tmp_2[i][j] = tmp_2[i][j-1]+dmci[i][j]/s0c[j];
				else tmp_2[i][j] = dmci[i][j]/s0c[j];
				
			}
		}
		else {
			for(int i=0; i<n; i++) {
				if(j > 0) tmp_2[i][j] = tmp_2[i][j-1];
				else tmp_2[i][j] = 0.;
			}
		}
	}//endj
	
	for(int i=0; i<n; i++) {		
		for(int j=0; j<njp; j++) {
			for(int k=0; k<ncovc; k++) ht[k][0] = h[i][j][k];
			jnktmp = matmult(transpose(ht),ric0inv);		
			for(int k=0; k<ncovc; k++) jnk0[i][j][k] = jnktmp[0][k];		
		}
	}//end i
	
	for(int i=0; i<n; i++) {
		for(int j=0; j<n; j++) {
			for(int k=0; k<ncov; k++) temp[j][k] = 0.;
			tj = 0;
			for(int j2=0; j2<njp; j2++) {
				if(ts[j] == tjp[j2]) {
					tj = j2;
					break;
				}
			}//end j2
			for(int j2=0; j2<njp; j2++) {
				if(tj < j2) {
					jnk1 = 0.;
					for(int k=0; k<ncovc; k++) jnk1 += jnk0[j][tj][k]*tmp_1[i][k];
					Wc1 = jnk1+tmp_2[i][tj]*bzic[j];
					
					jnk1 = 0;
					for(int k=0; k<ncovc; k++) jnk1 += jnk0[j][j2][k]*tmp_1[i][k];	
					Wc2 = jnk1+tmp_2[i][j2]*bzic[j];

					for(int k=0; k<ncov; k++) temp[j][k] += (zs[j][k]-s1s0[j2][k])*wdm1[j][j2]*(Wc1-Wc2);
				}
			}//end j2
			for(int k=0; k<ncov; k++) psihat[i][k] += temp[j][k];
		}
	}//endi

	for(int i=0; i<n; i++) {
		for(int k=0; k<ncov; k++) {
			tempsum2[i][k] = etahat[i][k]+psihat[i][k];
 			for(int l=0; l<ncov; l++) {
				temp2[i][k][l] += (etahat[i][k]+psihat[i][k])*(etahat[i][l]+psihat[i][l]);
				omegasum[k][l] += omegahat[i][k][l];
			}
		} 	
	}//endi
	
	for(int i=0; i<n; i++) {
		for(int k=0; k<ncov; k++) {
			for(int l=0; l<ncov; l++) tempsum[k][l] += temp2[i][k][l];
		}
	}//endi
	
	det = Determinant(omegasum,ncov);
	for(int k=0; k<ncov; k++) {
		for(int l=0; l<ncov; l++) {
			omega0[k][l] = omegasum[k][l];
			sigma[k][l] = tempsum[k][l];
		}
	}
	if(fabs(det) > 0) {
		omega0inv = invert(omega0);	
		final = matmult(matmult(omega0inv,sigma),omega0inv);
		for(int k=0; k<ncov; k++) {
			for(int l=0; l<ncov; l++) vcov[k][l] = final[k][l];
		}
		for(int k=0; k<ncov; k++) var[k] = vcov[k][k];		
	}
	else {
		Rprintf("Singular matrix!\n");
		return;
	}

	// baseline estimation
	for( int j=0; j<njp; j++) {
		event_or_not = 0.;
		for(int i=0; i<n; i++) {
			if(ics[i]==1 && wdn[i][j]>0) {
				event_or_not = 1.;
				break;
			}
		}
	
		if(event_or_not ==1)
		{
			for(int k=0; k<ncov; k++) {
				if(count4 > 0) Ewdn[count4][k] = Ewdn[count4-1][k];
				else Ewdn[count4][k] = 0.;
			}
			for(int i=0; i<n; i++) {
				for(int k=0; k<ncov; k++) {
					if(s0[j] > 0.) Ewdn[count4][k] += s1s0[j][k]*wdn[i][j]/s0[j];
				}
			}
			
			for(int k=0; k<ncov; k++) EwdnI[k] = 0.;
			for(int k=0; k<ncov; k++) {
				for(int l=0; l<ncov; l++) EwdnI[k] += Ewdn[count4][l]*omega0inv[k][l];  	
			}//end k,l

			for(int i=0; i<n; i++) {	
				jnk4 = 0.;
				for(int k=0; k<ncov; k++) jnk4 += EwdnI[k]*tempsum2[i][k];
				W_lambda2[i][count4] = jnk4;
			}//end i
			count4++;			
		}	
	}//endj

	for(int i=0; i<n; i++) {
		for(int j=0; j<njp; j++) W_lambda[i][j] = W_lambda1[i][j]-W_lambda2[i][j];	
	}
	for(int i=0; i<n; i++) {
		for(int j=0; j<count4; j++) Wlambda[i*count4+j] = W_lambda[i][j];
	}
	for(int i=0; i<n; i++) {
		for(int k=0; k<ncov; k++) {
			jnk6 = 0.;
			for (int l=0; l<ncov; l++) jnk6 += omega0inv[k][l]*tempsum2[i][l];
			Wbeta[i*ncov+k] = jnk6;
		}		
	}
	for(int j=0; j<count4; j++) {
		jnk5 = 0.;
		for (int i=0; i<n; i++) jnk5 += pow(W_lambda[i][j],2); 
		lambda10sd[j] = sqrt(jnk5);
	}	
	
	//memory cleanup
	for(int i=0; i<n; i++) delete[] Wc[i];
	delete[] Wc;
	for(int i=0; i<n; i++) delete[] tmp_2[i];
	delete[] tmp_2;	
	for(int i=0; i<n; i++) delete[] W_lambda[i];
	delete[] W_lambda;		
	for(int i=0; i<n; i++) delete[] W_lambda1[i];
	delete[] W_lambda1;	
	for(int i=0; i<n; i++) delete[] W_lambda2[i];
	delete[] W_lambda2;	
	for(int i=0; i<n; i++) delete[] wdm1[i];
	delete[] wdm1;
	for(int i=0; i<n; i++) delete[] yci[i];
	delete[] yci;
	for(int i=0; i<n; i++) delete[] dmci[i];
	delete[] dmci;
	for(int i=0; i<ncov; i++) delete[] omegasum[i];
	delete[] omegasum;
	
	for(int i=0; i<n; i++) {
		for(int j=0; j<njp; j++) delete[] h[i][j];
		delete[] h[i];
	}
	delete[] h;
	for(int i=0; i<n; i++) {
		for(int j=0; j<njp; j++) delete[] jnk0[i][j];
		delete[] jnk0[i];
	}
	delete[] jnk0;
}


// Code for Stratified KM implementation
FG_KM_Strata::FG_KM_Strata(int N, int NCOV, int NSTRATA, float *TS, int *ICS, float *ZS, int *STRATA_VAR) {
	n = N;
	ncov = NCOV;
	nstrata = NSTRATA;
	ts = Array1D<float> (n,0.);
	ics = Array1D<int> (n,0);
	zs = Array2D<float> (n,ncov,0.);
	strata_var = Array1D<int> (n,0);
	
	for(int i=0; i<n; i++) {
		ts[i] = TS[i];
		ics[i] = ICS[i];
		strata_var[i] = STRATA_VAR[i];
	}
	for(int i=0; i<n; i++) 
		for(int j=0; j<ncov; j++) zs[i][j] = ZS[i*ncov+j];
		
//	for(int i=0; i<n; i++) Rprintf("%f %d\n",zs[i][0], strata_var[i]);
}

void FG_KM_Strata::riskst(float **wdn1, float **wy1, float **dnc, float **yc, float **dnci, int *njp, vector<float> &tjp) const{
	//variable and arrays definition
	vector<int> idxt(n);
	Array2D<float> gkm(n,nstrata,0.);

	//fit censoring distribution
	tjp[0] = ts[0];
	*njp += 1;
	idxt[0] = *njp;

	for(int i=1; i<n; i++) {
		if (ts[i] > tjp[*njp-1]) {
			tjp[*njp] = ts[i];
			*njp += 1;
		}
		idxt[i] = *njp-1;
	}
	
	//calculate weights based on risk set	
	for(int j=0; j<*njp; j++) {
		for(int m=0; m<nstrata; m++) {
			dnc[j][m] = 0.;
			yc[j][m] = 0.;
		}
		for(int i=0; i<n; i++) {
			dnci[i][j] = 0.;
          	if (ts[i] == tjp[j] && ics[i] == 0) dnci[i][j] = 1.;
			for(int m=0; m<nstrata; m++) {
				if (ts[i] >= tjp[j] && strata_var[i] == m) yc[j][m] += 1;	
				if (ts[i] == tjp[j] && ics[i] == 0 && strata_var[i] == m) dnc[j][m] += 1;
			}
		}
		for(int m=0; m<nstrata; m++) {
			if (fabs(yc[j][m]) > 0 && j == 0) gkm[j][m] = 1.-dnc[j][m]/yc[j][m];
			else if(abs(yc[j][m] > 0) && j > 0) gkm[j][m] = gkm[j-1][m]*(1.-dnc[j][m]/yc[j][m]);
			else gkm[j][m] = gkm[j-1][m]; //endif
		}
		
		for(int i=0; i<n; i++) {
			wdn1[i][j] = 0.;
			wy1[i][j] = 0.;
			
			if (ts[i] >= tjp[j]) {
				if (gkm[j][strata_var[i]] > 0) wy1[i][j] = 1.;
				if (ts[i] == tjp[j] && ics[i] == 1) {
					if (gkm[j][strata_var[i]] > 0) wdn1[i][j] = 1.;
				}//endif
			}//endif
			
			if (ts[i] < tjp[j] && ics[i] > 1) {
				if (gkm[idxt[i]][strata_var[i]] > 0) wy1[i][j] = gkm[j][strata_var[i]]/gkm[idxt[i]][strata_var[i]];
			}//endif	
		}//endi			
	}//end j
	
	return;
}	

void FG_KM_Strata::fgest1(int njp, float **wy, float **wdn, vector<float> &gfg1, vector<float> &gfg2, float *iflag, int nvar) const {
	//variable declaration
	float s0, btmp, det;

	vector<float> bzi(n);
	vector<float> s1(nvar);
	vector<float> s1s0(nvar);

	Array2D<float> s2(nvar,nvar,0.);
	Array2D<float> s2s0(nvar,nvar,0.);
	Array2D<float> s1s0sq(nvar,nvar,0.);
		
	Array2D<float> ubeta(nvar,1,0.);
	Array2D<float> rbeta0(nvar,nvar,0.);
	Array2D<float> rbeta0inv(nvar,nvar,0.);
	Array2D<float> diff(nvar,1,0.);
	
	float **rbeta = new float*[nvar];
	for(int i=0; i<nvar; i++) rbeta[i] = new float[nvar];	
	
	//initialization
	for(int k=0; k<nvar; k++) {
		ubeta[k][0] = 0.;
		s1s0[k] = 0.;
		for(int l=0; l<nvar; l++) {
			rbeta[k][l] = 0.;
			s2s0[k][l] = 0.;
			s1s0sq[k][l] = 0.;
		}
	}//endk

	//calculation
	for (int i=0; i<n; i++) {
		btmp = 0.;	
		for(int k=0; k<nvar; k++) btmp += gfg1[k]*zs[i][k];		
		bzi[i] = exp(btmp);
	}//endi
	
	for(int j=0; j<njp; j++) {
		s0 = 0.;
		for(int k=0; k<nvar; k++) {
			s1[k] = 0.;
			for(int l=0; l<nvar; l++) {
				s2[k][l] = 0.;			
			}
		}//end k

		for (int i=0; i<n; i++) {
			s0 += wy[i][j]*bzi[i];
			for(int k=0; k<nvar; k++) s1[k] += wy[i][j]*zs[i][k]*bzi[i];
			for(int k=0; k<nvar; k++) {
				for(int l=0; l<nvar; l++) s2[k][l] += wy[i][j]*zs[i][k]*zs[i][l]*bzi[i];
			}
		}//end i

		if (fabs(s0) > 0) {
			for(int k=0; k<nvar; k++) s1s0[k] = s1[k]/s0;
          	for(int k=0; k<nvar; k++) {
				for(int l=0; l<nvar; l++) s2s0[k][l] = s2[k][l]/s0;
			}
		}
        else {
			for(int k=0; k<nvar; k++) s1s0[k] = 0.;
          	for(int k=0; k<nvar; k++) {
				for(int l=0; l<nvar; l++) s2s0[k][l] = 0.;
			}
		}//endif

		for(int k=0; k<nvar; k++) {
			for(int l=0; l<nvar; l++) s1s0sq[k][l] = s1s0[k]*s1s0[l];
		}
		
		for(int i=0; i<n; i++) {
			for(int k=0; k<nvar; k++) ubeta[k][0] += (zs[i][k]-s1s0[k])*wdn[i][j]; 
			for(int k=0; k<nvar; k++) {
				for(int l=0; l<nvar; l++) rbeta[k][l] += (s2s0[k][l]-s1s0sq[k][l])*wdn[i][j];
			}
		}//end i
	}//end j	
	
	
	for(int k=0; k<nvar; k++)
		for(int l=0; l<nvar; l++)
			rbeta0[k][l] = rbeta[k][l];

	det = Determinant(rbeta,nvar);
	if(fabs(det) > 0) {
		rbeta0inv = invert(rbeta0);	
		diff =  matmult(rbeta0inv,ubeta);
		for(int k=0; k<nvar; k++) {
			gfg2[k] = gfg1[k]+diff[k][0];
		}
	}
	else {
		*iflag = 1.;
	}			
	
	//memory cleanup
	for(int i=0; i<nvar; i++) delete[] rbeta[i];
	delete[] rbeta;
	
	return;
}

void FG_KM_Strata::fgest2km(vector<float> &tjp, int njp, float **wy, float **wdn, float gfg1[], float **dnc, float **yc, float **dnci, 
					 float var[], int nvar, float lambda10[], float tbase[], float lambda10sd[], float Wlambda[], float Wbeta[]) const {
	//variable declaration	
	int count = 0, count2 = 0, tj = 0, count4 = 0;
	float wdnsum, btmp, event_or_not, det, jnk4, jnk5, jnk6;
	
	vector<float> bzi(n);
	vector<float> a10hat(njp);
	vector<float> s0(njp);
	vector<float> s1(nvar);
	vector<float> EwdnI(nvar);
	
	Array2D<float> achat(n,nstrata,0.);
	Array2D<float> etahat(n,nvar,0.);
	Array2D<float> psihat(n,nvar,0.);
	Array2D<float> s1s0(n,nvar,0.);
	Array2D<float> temp(n,nvar,0.);
	Array2D<float> tempsum2(n,nvar,0.);
	Array2D<float> Ewdn(njp,nvar,0.);
	
	Array2D<float> omega0(nvar,nvar,0.);
	Array2D<float> omega0inv(nvar,nvar,0.);
	Array2D<float> sigma(nvar,nvar,0.);
	Array2D<float> final(nvar,nvar,0.);
	
	Array2D<float> s2(nvar,nvar,0.);
	Array2D<float> s2s0(nvar,nvar,0.);
	Array2D<float> s1s0sq(nvar,nvar,0.);
	Array2D<float> vcov(nvar,nvar,0.);
	Array2D<float> tempsum(nvar,nvar,0.);
	
	Array3D<float> temp2(n,nvar,nvar,0.);
	Array3D<float> omegahat(n,nvar,nvar,0.);
	
	float **wdm1 = new float*[n];
	for(int i=0; i<n; i++) wdm1[i] = new float[njp];	
	float **yci = new float*[n];
	for(int i=0; i<n; i++) yci[i] = new float[njp];		
	float **W_lambda = new float*[n];
	for(int i=0; i<n; i++) W_lambda[i] = new float[njp];	
	float **W_lambda1 = new float*[n];
	for(int i=0; i<n; i++) W_lambda1[i] = new float[njp];
	float **W_lambda2 = new float*[n];
	for(int i=0; i<n; i++) W_lambda2[i] = new float[njp];
	float **omegasum = new float*[nvar];
	for(int i=0; i<nvar; i++) omegasum[i] = new float[nvar];

	float ***dmci = new float**[n];
	for(int i=0; i<n; i++) dmci[i] = new float*[njp];
	for(int i=0; i<n; i++) {
		for(int j=0; j<njp; j++) dmci[i][j] = new float[nstrata];
	}

	float ***Wc = new float**[n];
	for(int i=0; i<n; i++) Wc[i] = new float*[njp];
	for(int i=0; i<n; i++) {
		for(int j=0; j<njp; j++) Wc[i][j] = new float[nstrata];
	}
	
	//initializaztion
	for(int i=0; i<n; i++) {
		btmp = 0.;	
		for(int k=0; k<nvar; k++) btmp += gfg1[k]*zs[i][k];		
		bzi[i] = exp(btmp);
		for(int k=0; k<nvar; k++) {
			etahat[i][k] = 0.;
			psihat[i][k] = 0.;
			for(int l=0; l<nvar; l++) omegahat[i][k][l]=0.;
		}
	}	
	
	for(int i=0; i<n; i++) {
		for(int k=0; k<nvar; k++) {
			for(int l=0; l<nvar; l++) temp2[i][k][l] = 0.;	
		}
	}
	
	for(int k=0; k<nvar; k++) 
		for(int l=0; l<nvar; l++) omegasum[k][l] = 0.;	

	for(int i=0; i<n; i++) {
		for(int k=0; k<nvar; k++) tempsum2[i][k] = 0.;
	}		
	
	for(int i=0; i<n; i++) {
		for(int k=0; k<nvar; k++) {
			for(int l=0; l<nvar; l++) tempsum[k][l] = 0.;
		}
	}
	
	for(int k=0; k<nvar; k++) EwdnI[k] = 0.;
	
	//calculation
	for(int j=0; j<njp; j++) {
		s0[j] = 0.;
		for(int k=0; k<nvar; k++) {
			s1[k] = 0.;
			for(int l=0; l<nvar; l++) {
				s2[k][l] = 0.;			
			}
		}//end k,l
		wdnsum = 0.;	
		for(int i=0; i<n; i++) {
			s0[j] += wy[i][j]*bzi[i];
			for(int k=0; k<nvar; k++) s1[k] += wy[i][j]*zs[i][k]*bzi[i];
			for(int k=0; k<nvar; k++) {
				for(int l=0; l<nvar; l++) s2[k][l] += wy[i][j]*zs[i][k]*zs[i][l]*bzi[i];
			}
			wdnsum += wdn[i][j];
		}//endi
		
		if(fabs(s0[j]) > 0) {	
			for(int k=0; k<nvar; k++) s1s0[j][k] = s1[k]/s0[j];
          	for(int k=0; k<nvar; k++) {
				for(int l=0; l<nvar; l++) s2s0[k][l] = s2[k][l]/s0[j];
			}
			if(j > 0) a10hat[j] = a10hat[j-1]+wdnsum/s0[j];
			else a10hat[j] = wdnsum/s0[j];
		}
		else {
			for(int k=0; k<nvar; k++) s1s0[j][k] = 0.;
          	for(int k=0; k<nvar; k++) {
				for(int l=0; l<nvar; l++) s2s0[k][l] = 0.;
			}	
			if(j > 0) a10hat[j] = a10hat[j-1];
			else a10hat[j] = 0.;
		}//endif
	
		//estimation of baseline hazard
		event_or_not = 0.;
		for(int i=0; i<n; i++) {
			if(ics[i]==1 && wdn[i][j]>0) {
				event_or_not = 1.;
				break;
			}
		}
		if(event_or_not == 1.) {
			tbase[count] = tjp[j];
			if(count > 0 && s0[j] > 0) lambda10[count] = lambda10[count-1]+wdnsum/s0[j];
			else if(s0[j] > 0) lambda10[count] = wdnsum/s0[j];
			count++;
		}//endif
		
		for(int k=0; k<nvar; k++) {
			for(int l=0; l<nvar; l++) s1s0sq[k][l] = s1s0[j][k]*s1s0[j][l];
		}
		
		for(int i=0; i<n; i++) {
			if(j > 0) wdm1[i][j] = wdn[i][j]-wy[i][j]*bzi[i]*(a10hat[j]-a10hat[j-1]);
			else wdm1[i][j] = wdn[i][j]-wy[i][j]*bzi[i]*a10hat[j];
			for(int k=0; k<nvar; k++) etahat[i][k] += (zs[i][k]-s1s0[j][k])*wdm1[i][j];
			for(int k=0; k<nvar; k++) {
				for(int l=0; l<nvar; l++) omegahat[i][k][l] += (s2s0[k][l]-s1s0sq[k][l])*wdn[i][j];
			}//end k,l
		}//endi			
		
		if(event_or_not == 1.) {         //calculate W_lambda1
			for(int i=0; i<n; i++) {
				if(count2 > 0 && s0[j] > 0) W_lambda1[i][count2] = W_lambda1[i][count2-1]+wdm1[i][j]/s0[j];
				else if( count2 == 0 && s0[j] > 0) W_lambda1[i][count2] = wdm1[i][j]/s0[j];
			}
			count2++;
		}	
	}//endj
	
	//calculate dmci & yci
	for(int j=0; j<njp; j++) {
		for(int m=0; m<nstrata; m++) {
			if(j > 0 && yc[j][m] > 0) achat[j][m] = achat[j-1][m]+dnc[j][m]/yc[j][m];
			else if(j > 0 && yc[j][m] == 0) achat[j][m] = achat[j-1][m];
			else if(j ==0 && yc[j][m] > 0) achat[j][m] = dnc[j][m]/yc[j][m];
		}
		for(int i=0; i<n; i++) {
			yci[i][j] = 0.;
			if (ts[i] >= tjp[j]) yci[i][j] = 1.;
			for(int m=0; m<nstrata; m++) {
				if(j > 0 && strata_var[i] == m) dmci[i][j][m] = dnci[i][j]-yci[i][j]*(achat[j][m]-achat[j-1][m]);   
				else dmci[i][j][m] = 0.;
			}
		}
	}//endj
	
	//calculate Wc matrix and variance
	for(int j=0; j<njp; j++) {
		for(int i=0; i<n; i++) {
			for(int m=0; m<nstrata; m++) {
				if(j > 0 && yc[j][m] > 0) Wc[i][j][m] = Wc[i][j-1][m]+dmci[i][j][m]/yc[j][m];
				else if(j > 0 && yc[j][m] == 0) Wc[i][j][m] = Wc[i][j-1][m];
				else if(j == 0 && yc[j][m] > 0) Wc[i][j][m] = dmci[i][j][m]/yc[j][m];
			}
		}
	}//endj	
	for(int i=0; i<n; i++) {
		for(int j=0; j<n; j++) {
			for(int k=0; k<nvar; k++) temp[j][k] = 0.;
			for(int j2=0; j2<njp; j2++) {
				if(ts[j] == tjp[j2]) {
					tj = j2;
					break;
				}
			}//end j2
			for(int j2=0; j2<njp; j2++) {
				if(tj < j2) {
					for(int k=0; k<nvar; k++) temp[j][k] += (zs[j][k]-s1s0[j2][k])*wdm1[j][j2]*(Wc[i][tj][strata_var[i]]-Wc[i][j2][strata_var[i]]);
				}
			}//end j2
			for(int k=0; k<nvar; k++) psihat[i][k] += temp[j][k];
		}
	}//end i,j

	for(int i=0; i<n; i++) {
		for(int k=0; k<nvar; k++) {
			tempsum2[i][k] = etahat[i][k]+psihat[i][k];                      
 			for(int l=0; l<nvar; l++) {
				temp2[i][k][l] += (etahat[i][k]+psihat[i][k])*(etahat[i][l]+psihat[i][l]);
				omegasum[k][l] += omegahat[i][k][l];
			}
		} 	
	}//endi	
	
	for(int i=0; i<n; i++) {
		for(int k=0; k<nvar; k++) {
			for(int l=0; l<nvar; l++) tempsum[k][l] += temp2[i][k][l];
		}
	}	

	det = Determinant(omegasum,nvar);
	for(int k=0; k<nvar; k++) {
		for(int l=0; l<nvar; l++) {
			omega0[k][l] = omegasum[k][l];
			sigma[k][l] = tempsum[k][l];
			if(isnan(sigma[k][l])) {
				Rprintf("C++ inner mistake.\n");
				return;
			}
		}
	}
	
	//cout << sigma << endl;			
	if(fabs(det) > 0.0001) {
		omega0inv = invert(omega0);		
		final = matmult(matmult(omega0inv,sigma),omega0inv);
		for(int k=0; k<nvar; k++) {
			for(int l=0; l<nvar; l++) vcov[k][l] = final[k][l];
		}
		for(int k=0; k<nvar; k++) var[k] = vcov[k][k];		
	}
	else {
		Rprintf("Singular matrix!\n");
		return;
	}
	
	//variance for baseline hazard
	for(int j=0; j<njp; j++) {
		event_or_not = 0.;
		for(int i=0; i<n; i++) {
			if(ics[i]==1 && wdn[i][j]>0) {
				event_or_not = 1.;
				break;
			}
		}
		if(event_or_not ==1.) {
			for(int k=0; k<nvar; k++) {
				if(count4 > 0) Ewdn[count4][k] = Ewdn[count4-1][k];
				else if(count4 == 0) Ewdn[count4][k] = 0.;
			}
			for(int i=0; i<n; i++) {
				if(s0[j] > 0) {
					for(int k=0; k<nvar; k++) Ewdn[count4][k] += s1s0[j][k]*wdn[i][j]/s0[j];
				}
			}
			
			for(int k=0; k<nvar; k++) EwdnI[k] = 0.;
			for(int k=0; k<nvar; k++) {
				for(int l=0; l<nvar; l++) EwdnI[k] += Ewdn[count4][l]*omega0inv[l][k];
			}

			for(int i=0; i<n; i++) {	
				jnk4 = 0.;
				for(int k=0; k<nvar; k++) jnk4 += EwdnI[k]*tempsum2[i][k];
				W_lambda2[i][count4] = jnk4;
			}	
			count4++;			
		}
	}//endj	
	
	for(int i=0; i<n; i++) {
		for(int j=0; j<count4; j++) W_lambda[i][j] = W_lambda1[i][j]-W_lambda2[i][j];
	}
	for(int i=0; i<n; i++) {
		for(int j=0; j<count4; j++) Wlambda[i*count4+j] = W_lambda[i][j];
	}	
	
	for(int i=0; i<n; i++) {
		for(int k=0; k<nvar; k++) {
			jnk6 = 0.;
			for (int l=0; l<nvar; l++) jnk6 += omega0inv[k][l]*tempsum2[i][l];
			Wbeta[i*nvar+k] = jnk6;
		}		
	}
	for(int j=0; j<count4; j++) {
		jnk5 = 0.;
		for(int i=0; i<n; i++) jnk5 += pow((double)W_lambda[i][j],2.);
		lambda10sd[j] = sqrt(jnk5);
	}	
	
	//memory cleanup
	for(int i=0; i<n; i++) delete[] wdm1[i];
	delete[] wdm1;
	for(int i=0; i<n; i++) delete[] yci[i];
	delete[] yci;
	for(int i=0; i<n; i++) delete[] W_lambda[i];
	delete[] W_lambda;
	for(int i=0; i<n; i++) delete[] W_lambda1[i];
	delete[] W_lambda1;
	for(int i=0; i<n; i++) delete[] W_lambda2[i];
	delete[] W_lambda2;
	for(int i=0; i<nvar; i++) delete[] omegasum[i];
	delete[] omegasum;
	
	for(int i=0; i<n; i++) {
		for(int j=0; j<njp; j++) delete[] dmci[i][j];
		delete[] dmci[i];
	}
	delete[] dmci;

	for(int i=0; i<n; i++) {
		for(int j=0; j<njp; j++) delete[] Wc[i][j];
		delete[] Wc[i];
	}
	delete[] Wc;
}


// Code for Stratified Cox implementation
FG_COX_Strata::FG_COX_Strata(int N, int NCOV, int NCOVC, int NSTRATA, float *TS, int *ICS, float *ZS, float *ZCS, int *STRATA_VAR) {
	n = N;
	ncov = NCOV;
	ncovc = NCOVC;
	nstrata = NSTRATA;
	ts = Array1D<float> (n,0.);
	ics = Array1D<int> (n,0.);
	zs = Array2D<float> (n,ncov,0.);
	zcs = Array2D<float> (n,ncovc,0.);
	strata_var = Array1D<int> (n,0);
	
	for(int i=0; i<n; i++) {
		ts[i] = TS[i];
		ics[i] = ICS[i];
		strata_var[i] = STRATA_VAR[i];
	}
	for(int i=0; i<n; i++) {
		for(int j=0; j<ncov; j++) zs[i][j] = ZS[i*ncov+j];
		for(int j=0; j<ncovc; j++) zcs[i][j] = ZCS[i*ncovc+j];
	}
}

void FG_COX_Strata::riskst(float **wdn3, float **wy3, float *iflagc, vector<float> &betac, float **dnc, float **yc, float **dnci, int *njp, vector<float> &tjp, float **ric, float **rlamc0, float ***s1s0c, int nvar, float *censdet) const{
	//variable and arrays definition
	int nit = 0;
	float tmp, max_error=1.;
	
	vector<int> idxt(n);
	vector<float> s0(nstrata);
	vector<float> bzi(n);
	vector<float> error(nvar);
	vector<float> bc1(nvar);
	vector<float> bc2(nvar);
	
	Array2D<float> s1(nvar,nstrata,0.);
	
	Array3D<float> s2(nvar,nvar,nstrata);
	Array3D<float> s1s0sq(nvar,nvar,nstrata);
	Array3D<float> s2s0(nvar,nvar,nstrata);

	void ccox_strata(int n, const Array1D<float> &ts, const Array1D<int> &ics, const Array2D<float> &zs, vector<float> &tjp, int njp, vector<float> &bc1, vector<float> &bc2, float **dnci, float *iflagc, int nvar, float *censdet, int nstrata, const Array1D<int> &strata_var);
	
	//dynamic arrays allocation
	float **gsurv = new float*[n];
	for(int i=0; i<n; i++) gsurv[i] = new float[n];	

	//fit censoring distribution
	for(int k=0; k<nvar; k++) {
		bc1[k] = 0.;
		bc2[k] = 0.;
	}

	tjp[0] = ts[0];
	*njp += 1;
	idxt[0] = *njp;

	for(int i=1; i<n; i++) {
		if (ts[i] > tjp[*njp-1]) {
			tjp[*njp] = ts[i];
			*njp += 1;
		}
		idxt[i] = *njp-1;
	}
	
	while(max_error > 0.00001 && nit < 50) {
		ccox_strata(n, ts, ics, zcs, tjp, *njp, bc1, bc2, dnci, iflagc, nvar, censdet, nstrata, strata_var);
		nit++;
		if(*iflagc==1) return;
		else {
			max_error = fabs(bc2[0]-bc1[0]);
			for(int k=1; k<nvar; k++) {
				error[k] = fabs(bc2[k]-bc1[k]);		
				if(error[k] > max_error) max_error = error[k];
			}
			for(int k=0; k<nvar; k++) bc1[k] = bc2[k];		
		}//endif	
	}//end while

	if (nit==50) {
		*iflagc = 1;
		return;
	}	
	for(int k=0; k<nvar; k++) betac[k]=bc2[k];	

	//calculate weights based on risk set
	for (int i=0; i<n; i++) {
		tmp = 0.;	
		for(int k=0; k<nvar; k++) {
			tmp += betac[k]*zcs[i][k];
		}
		bzi[i] = exp(tmp);
	}
	for(int k=0; k<nvar; k++)
		for(int l=0; l<nvar; l++) ric[k][l] = 0.;	

	for(int j=0; j<*njp; j++) {
		for(int m=0; m<nstrata; m++) {
			dnc[j][m] = 0.;
			yc[j][m] = 0.;
			s0[m] = 0.;
			for(int k=0; k<nvar; k++) {
				s1[k][m] = 0.;
				for(int l=0; l<nvar; l++) {
					s2[k][l][m] = 0.;			
				}	
			}			
			
			for(int i=0; i<n; i++) {
				if (ts[i] >= tjp[j] && strata_var[i] == m) {
					yc[j][m] += 1;
					s0[m] += bzi[i];
					for(int k=0; k<nvar; k++) s1[k][m] += zcs[i][k]*bzi[i];
					for(int k=0; k<nvar; k++) {
						for(int l=0; l<nvar; l++) s2[k][l][m] += zcs[i][k]*zcs[i][l]*bzi[i]; 
					}
				}//endif	
				if (ts[i] == tjp[j] && ics[i] == 0 && strata_var[i] == m) dnc[j][m] += 1;
			}

			if (fabs(s0[m]) > 0.) {
				for(int k=0; k<nvar; k++) s1s0c[j][k][m] = s1[k][m]/s0[m];
				for(int k=0; k<nvar; k++) {
					for(int l=0; l<nvar; l++) s2s0[k][l][m] = s2[k][l][m]/s0[m];
				}
				if(j == 0) rlamc0[j][m] = dnc[j][m]/s0[m];
				else rlamc0[j][m] = rlamc0[j-1][m] + dnc[j][m]/s0[m];
			}			
			else {
				for(int k=0; k<nvar; k++) s1s0c[j][k][m] = 0.;
				for(int k=0; k<nvar; k++) {
					for(int l=0; l<nvar; l++) s2s0[k][l][m] = 0.;
				}
				if(j > 0) rlamc0[j][m] = rlamc0[j-1][m];
			} //endif

			for(int k=0; k<nvar; k++) {
				for(int l=0; l<nvar; l++) s1s0sq[k][l][m] = s1s0c[j][k][m]*s1s0c[j][l][m];
			}		
			for(int k=0; k<nvar; k++) {
				for(int l=0; l<nvar; l++) ric[k][l] += (s2s0[k][l][m] - s1s0sq[k][l][m])*dnc[j][m];
			}

			for(int i=0; i<n; i++) {
				if(strata_var[i] == m) {
					gsurv[i][j] = exp(0.-rlamc0[j][m]*bzi[i]);
					wdn3[i][j] = 0.;
					wy3[i][j] = 0.;
					if (ts[i] >= tjp[j]) {
						if (gsurv[i][j] > 0) wy3[i][j] = 1.;
						if (ts[i] == tjp[j] && ics[i] == 1) {
							if (gsurv[i][j] > 0) wdn3[i][j] = 1.;
						}
					}
					if (ts[i] < tjp[j] && ics[i] > 1) {
						if (gsurv[i][idxt[i]] > 0) wy3[i][j] = gsurv[i][j]/gsurv[i][idxt[i]];
					}
				}//endif
			}//endi			
		}//endm
	}//end j
	
	//memory cleanup
	for(int i=0; i<n; i++) delete[] gsurv[i];
	delete[] gsurv;		
}

void ccox_strata(int n, const Array1D<float> &ts, const Array1D<int> &ics, const Array2D<float> &zs, vector<float> &tjp, int njp, vector<float> &bc1, vector<float> &bc2, float **dnci, float *iflagc, int nvar, float *censdet, int nstrata, const Array1D<int> &strata_var) {
	//variable declaration
	float btmp, det;

	vector<float> s0(nstrata);
	
	Array2D<float> s1(nvar,nstrata,0.);
	Array2D<float> s1s0(nvar,nstrata,0.);
	
	Array2D<float> uc(nvar,1,0.);
	Array2D<float> rc0(nvar,nvar,0.);
	Array2D<float> rc0inv(nvar,nvar,0.);
	Array2D<float> diff(nvar,1,0.);
	
	Array3D<float> s2(nvar,nvar,nstrata,0.);
	Array3D<float> s2s0(nvar,nvar,nstrata,0.);
	Array3D<float> s1s0sq(nvar,nvar,nstrata,0.);
	
	float **rc = new float*[nvar];
	for(int i=0; i<nvar; i++) rc[i] = new float[nvar];
	
	//initialization
	for(int k=0; k<nvar; k++) {
		uc[k][0] = 0.;
		for(int m=0; m<nstrata; m++) s1s0[k][m] = 0.;
		for(int l=0; l<nvar; l++) {
			rc[k][l] = 0.;
			for(int m=0; m<nstrata; m++) {
				s2s0[k][l][m] = 0.;
				s1s0sq[k][l][m] = 0.;
			}
		}
	}			  
		  
	//calculation
	for (int j=0; j<njp; j++) {
		for(int m=0; m<nstrata; m++) {
			s0[m] = 0.;
			for(int k=0; k<nvar; k++) {
				s1[k][m] = 0.;
				for(int l=0; l<nvar; l++) {
					s2[k][l][m] = 0.;			
				}
			}//end k
			
			for(int i=0; i<n; i++) {
				if(strata_var[i] == m) {
					btmp = 0.;
					double tmp = 0.;
					for(int k=0; k<nvar; k++) tmp += bc1[k]*zs[i][k];
					btmp = exp(tmp);
						
					if (ts[i] >= tjp[j]) {
						s0[m] += btmp;
						for(int k=0; k<nvar; k++) s1[k][m] += zs[i][k]*btmp;
						for(int k=0; k<nvar; k++) {
							for(int l=0; l<nvar; l++) s2[k][l][m] += zs[i][k]*zs[i][l]*btmp;
						}
					}//endif				
				}
			}//endi
			
			if(fabs(s0[m]) > 0.) {
				for(int k=0; k<nvar; k++) s1s0[k][m] = s1[k][m]/s0[m];
				for(int k=0; k<nvar; k++)
					for(int l=0; l<nvar; l++)
						s2s0[k][l][m] = s2[k][l][m]/s0[m];
			}
			else {
				for(int k=0; k<nvar; k++) s1s0[k][m] = 0.;
				for(int k=0; k<nvar; k++)
					for(int l=0; l<nvar; l++)
							s2s0[k][l][m] = 0.;
			}//endif

			for(int k=0; k<nvar; k++) {
				for(int l=0; l<nvar; l++) s1s0sq[k][l][m] = s1s0[k][m]*s1s0[l][m];
			}			
		
			for(int i=0; i<n; i++) {
				if(strata_var[i] == m) {
					dnci[i][j] = 0.;
					if (ts[i] == tjp[j] && ics[i] == 0) dnci[i][j] = 1.;
					for(int k=0; k<nvar; k++) {
						uc[k][0] += (zs[i][k] - s1s0[k][m])*dnci[i][j];	
					}
					for(int k=0; k<nvar; k++) {
						for(int l=0; l<nvar; l++) {
							rc[k][l] += (s2s0[k][l][m] - s1s0sq[k][l][m])*dnci[i][j];
						}
					}
				}	
			}//endi
		}//endm
	}//end j			  

	for(int k=0; k<nvar; k++)
		for(int l=0; l<nvar; l++)
			rc0[k][l] = rc[k][l];	

	det = Determinant(rc,nvar);
	
	if(fabs(det) > 0) {
		if(fabs(det) < 0.00001) {
			censdet[0] = 0.;
			iflagc[0] = 1.;
			return;
		}
		rc0inv = invert(rc0);	
		diff =  matmult(rc0inv,uc);
		for(int k=0; k<nvar; k++) {
			bc2[k] = bc1[k]+diff[k][0];
		}
	}
	else {
		*iflagc = 1.;
		return;
	}
	
	//memory cleanup
	for(int i=0; i<nvar; i++) delete[] rc[i];
	delete[] rc;
			  
}

void FG_COX_Strata::fgest1(int njp, float **wy, float **wdn, vector<float> &gfg1, vector<float> &gfg2, float *iflag, int nvar) const {
	//variable declaration
	float s0, btmp, det;
	
	vector<float> bzi(n);
	vector<float> s1(nvar);
	vector<float> s1s0(nvar);
	
	Array2D<float> s2(nvar,nvar,0.);
	Array2D<float> s2s0(nvar,nvar,0.);
	Array2D<float> s1s0sq(nvar,nvar,0.);
	
	Array2D<float> ubeta(nvar,1,0.);
	Array2D<float> rbeta0(nvar,nvar,0.);
	Array2D<float> rbeta0inv(nvar,nvar,0.);
	Array2D<float> diff(nvar,1,0.);
	
	float **rbeta = new float*[nvar];
	for(int i=0; i<nvar; i++) rbeta[i] = new float[nvar];	
	
	//initialization
	for(int k=0; k<nvar; k++) {
		ubeta[k][0] = 0.;
		s1s0[k] = 0.;
		for(int l=0; l<nvar; l++) {
			rbeta[k][l] = 0.;
			s2s0[k][l] = 0.;
			s1s0sq[k][l] = 0.;
		}
	}//endk

	//calculation
	for (int i=0; i<n; i++) {
		btmp = 0.;	
		for(int k=0; k<nvar; k++) btmp += gfg1[k]*zs[i][k];		
		bzi[i] = exp(btmp);
	}//endi
	
	for(int j=0; j<njp; j++) {
		s0 = 0.;
		for(int k=0; k<nvar; k++) {
			s1[k] = 0.;
			for(int l=0; l<nvar; l++) {
				s2[k][l] = 0.;			
			}
		}//end k

		for (int i=0; i<n; i++) {
			s0 += wy[i][j]*bzi[i];
			for(int k=0; k<nvar; k++) s1[k] += wy[i][j]*zs[i][k]*bzi[i];
			for(int k=0; k<nvar; k++) {
				for(int l=0; l<nvar; l++) s2[k][l] += wy[i][j]*zs[i][k]*zs[i][l]*bzi[i];
			}
		}//end i

		if (fabs(s0) > 0) {
			for(int k=0; k<nvar; k++) s1s0[k] = s1[k]/s0;
          	for(int k=0; k<nvar; k++) {
				for(int l=0; l<nvar; l++) s2s0[k][l] = s2[k][l]/s0;
			}
		}
        else {
			for(int k=0; k<nvar; k++) s1s0[k] = 0.;
          	for(int k=0; k<nvar; k++) {
				for(int l=0; l<nvar; l++) s2s0[k][l] = 0.;
			}
		}//endif

		for(int k=0; k<nvar; k++) {
			for(int l=0; l<nvar; l++) s1s0sq[k][l] = s1s0[k]*s1s0[l];
		}
		
		for(int i=0; i<n; i++) {
			for(int k=0; k<nvar; k++) ubeta[k][0] += (zs[i][k]-s1s0[k])*wdn[i][j]; 
			for(int k=0; k<nvar; k++) {
				for(int l=0; l<nvar; l++) rbeta[k][l] += (s2s0[k][l]-s1s0sq[k][l])*wdn[i][j];
			}
		}//end i
	}//end j	

	for(int k=0; k<nvar; k++)
		for(int l=0; l<nvar; l++)
			rbeta0[k][l] = rbeta[k][l];

	det = Determinant(rbeta,nvar);
	if(fabs(det) > 0) {
		rbeta0inv = invert(rbeta0);	
		diff =  matmult(rbeta0inv,ubeta);
		for(int k=0; k<nvar; k++) {
			gfg2[k] = gfg1[k]+diff[k][0];
		}
	}
	else {
		*iflag = 1.;
	}			
	
	//memory cleanup
	for(int i=0; i<nvar; i++) delete[] rbeta[i];
	delete[] rbeta;
	
	return;
}

void FG_COX_Strata::fgest2cox(vector<float> &tjp, int njp, vector<float> &betac, float **ric, float **rlamc0, float ***s1s0c, float **wy, float **wdn, float gfg1[], float **dnc, float **yc, float **dnci, float var[], float lambda10[], float tbase[], float lambda10sd[], float Wlambda[], float Wbeta[], float *iflagric) const {
	//variable declaration
	int count = 0, count2 = 0, count4 = 0, tj;
	float wdnsum, btmp, btmp2, event_or_not, det, Wc1 = 0., Wc2 = 0., jnk1 = 0., jnk4, jnk5, jnk6;
	
	vector<float> bzi(n);
	vector<float> bzic(n);
	vector<float> a10hat(n);
	vector<float> s0(njp);
	vector<float> s1(ncov);
	vector<float> EwdnI(ncov);
	
	Array2D<float> s0c(njp,nstrata,0.);
	Array2D<float> s1s0(njp,ncov,0.);
	Array2D<float> s2(ncov,ncov,0.);
	Array2D<float> s2s0(ncov,ncov,0.);
	Array2D<float> s1s0sq(ncov,ncov,0.);
	Array2D<float> temp(njp,ncov,0.);
	Array2D<float> Ewdn(njp,ncov,0.);
	
	Array2D<float> etahat(n,ncov,0.);
	Array2D<float> psihat(n,ncov,0.);
	Array2D<float> tmp_1(n,ncovc,0.);
	Array2D<float> tempsum(ncov,ncov,0.);
	Array2D<float> tempsum2(n,ncov,0.);
	Array2D<float> vcov(ncov,ncov,0.);
	
	Array2D<float> ric0(ncovc,ncovc,0.);
	Array2D<float> ric0inv(ncovc,ncovc,0.);
	Array2D<float> ht(ncovc,1,0.);
	Array2D<float> jnktmp(1,ncovc,0.);
	Array2D<float> omega0(ncov,ncov,0.);
	Array2D<float> omega0inv(ncov,ncov,0.);
	Array2D<float> sigma(ncov,ncov,0.);
	Array2D<float> final(ncov,ncov,0.);
	
	Array3D<float> temp2(n,ncov,ncov,0.);
	Array3D<float> omegahat(n,ncov,ncov,0.);
	
	float **Wc = new float*[n];
	for(int i=0; i<n; i++) Wc[i] = new float[njp];	
	float **tmp_2 = new float*[n];
	for(int i=0; i<n; i++) tmp_2[i] = new float[njp];
	float **W_lambda = new float*[n];
	for(int i=0; i<n; i++) W_lambda[i] = new float[njp];
	float **W_lambda1 = new float*[n];
	for(int i=0; i<n; i++) W_lambda1[i] = new float[njp];
	float **W_lambda2 = new float*[n];
	for(int i=0; i<n; i++) W_lambda2[i] = new float[njp];
	float **wdm1 = new float*[n];
	for(int i=0; i<n; i++) wdm1[i] = new float[njp];	
	float **yci = new float*[n];
	for(int i=0; i<n; i++) yci[i] = new float[njp];	
	float **dmci = new float*[n];
	for(int i=0; i<n; i++) dmci[i] = new float[njp];	
	float **omegasum = new float*[ncov];
	for(int i=0; i<ncov; i++) omegasum[i] = new float[ncov];
	
	float ***h = new float**[n];
	for(int i=0; i<n; i++) h[i] = new float*[njp];
	for(int i=0; i<n; i++) {
		for(int j=0; j<njp; j++) h[i][j] = new float[ncovc];
	}
	float ***jnk0 = new float**[n];
	for(int i=0; i<n; i++) jnk0[i] = new float*[njp];
	for(int i=0; i<n; i++) {
		for(int j=0; j<njp; j++) jnk0[i][j] = new float[ncovc];
	}

	
	//basic calculation and initialization
	for(int i=0; i<n; i++) {
		btmp = 0.;
		btmp2 = 0.;	
		for(int k=0; k<ncov; k++) btmp += gfg1[k]*zs[i][k];
		for(int k=0; k<ncovc; k++) btmp2 += betac[k]*zcs[i][k];

		bzi[i] = exp(btmp);
		bzic[i] = exp(btmp2);
		for(int k=0; k<ncov; k++) {
			etahat[i][k] = 0.;
			psihat[i][k] = 0.;
			for(int l=0; l<ncov; l++) omegahat[i][k][l]=0.;
		}
	}//endi	
	
	for(int i=0; i<n; i++) {
		for(int k=0; k<ncovc; k++) 
			tmp_1[i][k] = 0.;
	}	

	for(int i=0; i<n; i++) {
		for(int k=0; k<ncov; k++) {
			for(int l=0; l<ncov; l++) temp2[i][k][l] = 0.;
		}
	}

	for(int k=0; k<ncov; k++) 
		for(int l=0; l<ncov; l++) omegasum[k][l] = 0.;

	for(int i=0; i<n; i++)
		for(int k=0; k<ncov; k++) tempsum2[i][k] = 0.;

	for(int k=0; k<ncov; k++)
		for(int l=0; l<ncov; l++) tempsum[k][l] = 0.;
	
	for(int k=0; k<ncov; k++) EwdnI[k] = 0.;	
	
	//calculation for variance
	// Major term
	for(int j=0; j<njp; j++) {
		s0[j] = 0.;
		for(int k=0; k<ncov; k++) {
			s1[k] = 0.;
			for(int l=0; l<ncov; l++) s2[k][l] = 0.;			
		}
		wdnsum = 0.;

		for(int i=0; i<n; i++) {
			s0[j] += wy[i][j]*bzi[i];
			for(int k=0; k<ncov; k++) s1[k] += wy[i][j]*zs[i][k]*bzi[i];
			for(int k=0; k<ncov; k++) {
				for(int l=0; l<ncov; l++) s2[k][l] += wy[i][j]*zs[i][k]*zs[i][l]*bzi[i];
			}
			wdnsum += wdn[i][j];
		}//endi	
		
		if (fabs(s0[j]) > 0) {	
			for(int k=0; k<ncov; k++) s1s0[j][k] = s1[k]/s0[j];
          	for(int k=0; k<ncov; k++) {
				for(int l=0; l<ncov; l++) s2s0[k][l] = s2[k][l]/s0[j];
			}
			if(j > 0) a10hat[j] = a10hat[j-1]+wdnsum/s0[j];
			else a10hat[j] = wdnsum/s0[j];
		}
		else {
			for(int k=0; k<ncov; k++) s1s0[j][k] = 0.;
          	for(int k=0; k<ncov; k++) {
				for(int l=0; l<ncov; l++) s2s0[k][l] = 0.;
			}	
			if(j > 0) a10hat[j] = a10hat[j-1];
			else a10hat[j] = 0.;
		}//endif

		event_or_not = 0.;
		for(int i=0; i<n; i++) {
			if(ics[i]==1 && wdn[i][j]>0) {
				event_or_not = 1.;
				break;
			}
		}		
		if(event_or_not == 1) {
			tbase[count] = tjp[j];
			if(count > 0) lambda10[count] = lambda10[count-1]+wdnsum/s0[j];
			else lambda10[count] = wdnsum/s0[j];
			count++;
		}//endif

		for(int k=0; k<ncov; k++) {
			for(int l=0; l<ncov; l++) s1s0sq[k][l] = s1s0[j][k]*s1s0[j][l];
		}

		for(int i=0; i<n; i++) {
			if(j > 0) wdm1[i][j] = wdn[i][j]-wy[i][j]*bzi[i]*(a10hat[j]-a10hat[j-1]);
			else wdm1[i][j] = wdn[i][j]-wy[i][j]*bzi[i]*a10hat[j];
			for(int k=0; k<ncov; k++) etahat[i][k] += (zs[i][k]-s1s0[j][k])*wdm1[i][j];
			for(int k=0; k<ncov; k++) {
				for(int l=0; l<ncov; l++) omegahat[i][k][l] += (s2s0[k][l]-s1s0sq[k][l])*wdn[i][j];
			}
		}//endi	
		
		if(event_or_not==1.) {
			for(int i=0; i<n; i++) {
				if(count2 > 0) W_lambda1[i][count2] = W_lambda1[i][count2-1]+wdm1[i][j]/s0[j];
				else W_lambda1[i][count2] = wdm1[i][j]/s0[j];
			}
			count2++;
		}				
	}//endj

	// Minor term
	det = Determinant(ric,ncovc);
	for(int k=0; k<ncovc; k++) {
		for(int l=0; l<ncovc; l++) ric0[k][l] = ric[k][l];
	}
	
	if(fabs(det) > 0.0001) ric0inv = invert(ric0);	
	else {
		iflagric[0] = 1.;
		return;
	}
	
	for(int j=0; j<njp; j++) {
		for(int i=0; i<n; i++) {
			yci[i][j] = 0.;
			if (ts[i] >= tjp[j]) yci[i][j] = 1.;
			if(j > 0) dmci[i][j] = dnci[i][j]-yci[i][j]*bzic[i]*(rlamc0[j][strata_var[i]]-rlamc0[j-1][strata_var[i]]);
			else dmci[i][j] = dnci[i][j]-yci[i][j]*bzic[i]*rlamc0[j][strata_var[i]];
		}
	}//endj

	for(int m=0; m<nstrata; m++) {
		for(int j=0; j<njp; j++) {
			s0c[j][m] = 0.;	
			for(int i=0; i<n; i++) {
				if(ts[i] >= tjp[j] && strata_var[i] == m) s0c[j][m] += bzic[i];
			}
		}
	}//endm

	for(int j=0; j<njp; j++) {
		for(int i=0; i<n; i++) {
			for(int k=0; k<ncovc; k++) {
				if(j > 0) h[i][j][k] = h[i][j-1][k]+bzic[i]*(zcs[i][k]-s1s0c[j][k][strata_var[i]])*(rlamc0[j][strata_var[i]]-rlamc0[j-1][strata_var[i]]);
				else h[i][j][k] = bzic[i]*(zcs[i][k]-s1s0c[j][k][strata_var[i]])*rlamc0[j][strata_var[i]];
			}
		}
	}//end j

	for(int i=0; i<n; i++) {
		for(int j=0; j<njp; j++) {
			for(int k=0; k<ncovc; k++) tmp_1[i][k] += (zcs[i][k]-s1s0c[j][k][strata_var[i]])*dmci[i][j];
		}
	}//end i	
		
	for(int j=0; j<njp; j++) {
		for(int i=0; i<n; i++) {
			if(s0c[j][strata_var[i]] > 0.) {
				if(j > 0) tmp_2[i][j] = tmp_2[i][j-1]+dmci[i][j]/s0c[j][strata_var[i]];
				else tmp_2[i][j] = dmci[i][j]/s0c[j][strata_var[i]];
			}
			else {
				if(j > 0) tmp_2[i][j] = tmp_2[i][j-1];
				else tmp_2[i][j] = 0.;		
			}
		}
	}//endj		

	for(int i=0; i<n; i++) {
		for(int j=0; j<njp; j++) {
			for(int k=0; k<ncovc; k++) ht[k][0] = h[i][j][k];
			jnktmp = matmult(transpose(ht),ric0inv);		
			for(int k=0; k<ncovc; k++) jnk0[i][j][k] = jnktmp[0][k];		
		}
	}//end i						

	for(int i=0; i<n; i++) {
		for(int j=0; j<n; j++) {
			for(int k=0; k<ncov; k++) temp[j][k] = 0.;
			tj = 0;
			for(int j2=0; j2<njp; j2++) {
				if(ts[j] == tjp[j2]) {
					tj = j2;
					break;
				}
			}//end j2
			for(int j2=0; j2<njp; j2++) {
				if(tj < j2) {
					jnk1 = 0.;
					for(int k=0; k<ncovc; k++) jnk1 += jnk0[j][tj][k]*tmp_1[i][k];
					Wc1 = jnk1+tmp_2[i][tj]*bzic[j];
				
					jnk1 = 0;
					for(int k=0; k<ncovc; k++) jnk1 += jnk0[j][j2][k]*tmp_1[i][k];	
					Wc2 = jnk1+tmp_2[i][j2]*bzic[j];

					for(int k=0; k<ncov; k++) temp[j][k] += (zs[j][k]-s1s0[j2][k])*wdm1[j][j2]*(Wc1-Wc2);
				}
			}//end j2
			for(int k=0; k<ncov; k++) psihat[i][k] += temp[j][k];
		}
	}//endi		


	// Variance for betas
	for(int i=0; i<n; i++) {
		for(int k=0; k<ncov; k++) {
			tempsum2[i][k] = etahat[i][k]+psihat[i][k];
 			for(int l=0; l<ncov; l++) {
				temp2[i][k][l] += (etahat[i][k]+psihat[i][k])*(etahat[i][l]+psihat[i][l]);
				omegasum[k][l] += omegahat[i][k][l];
			}
		} 	
	}//endi
	
	for(int i=0; i<n; i++) {
		for(int k=0; k<ncov; k++) {
			for(int l=0; l<ncov; l++) tempsum[k][l] += temp2[i][k][l];
		}
	}//endi
	
	det = Determinant(omegasum,ncov);
	for(int k=0; k<ncov; k++) {
		for(int l=0; l<ncov; l++) {
			omega0[k][l] = omegasum[k][l];
			sigma[k][l] = tempsum[k][l];
		}
	}
	if(fabs(det) > 0) {
		omega0inv = invert(omega0);	
		final = matmult(matmult(omega0inv,sigma),omega0inv);
		for(int k=0; k<ncov; k++) {
			for(int l=0; l<ncov; l++) vcov[k][l] = final[k][l];
		}
		for(int k=0; k<ncov; k++) var[k] = vcov[k][k];		
	}
	else {
		Rprintf("Singular matrix!\n");
		return;
	}

	// baseline estimation & variance
	for( int j=0; j<njp; j++) {
		event_or_not = 0.;
		for(int i=0; i<n; i++) {
			if(ics[i]==1 && wdn[i][j]>0) {
				event_or_not = 1.;
				break;
			}
		}
	
		if(event_or_not ==1)
		{
			for(int k=0; k<ncov; k++) {
				if(count4 > 0) Ewdn[count4][k] = Ewdn[count4-1][k];
				else Ewdn[count4][k] = 0.;
			}
			for(int i=0; i<n; i++) {
				for(int k=0; k<ncov; k++) {
					if(s0[j] > 0.) Ewdn[count4][k] += s1s0[j][k]*wdn[i][j]/s0[j];
				}
			}
			
			for(int k=0; k<ncov; k++) EwdnI[k] = 0.;
			for(int k=0; k<ncov; k++) {
				for(int l=0; l<ncov; l++) EwdnI[k] += Ewdn[count4][l]*omega0inv[k][l];  	
			}//end k,l

			for(int i=0; i<n; i++) {	
				jnk4 = 0.;
				for(int k=0; k<ncov; k++) jnk4 += EwdnI[k]*tempsum2[i][k];
				W_lambda2[i][count4] = jnk4;
			}//end i
			count4++;			
		}	
	}//endj

	for(int i=0; i<n; i++) {
		for(int j=0; j<njp; j++) W_lambda[i][j] = W_lambda1[i][j]-W_lambda2[i][j];	
	}
	for(int i=0; i<n; i++) {
		for(int j=0; j<count4; j++) Wlambda[i*count4+j] = W_lambda[i][j];
	}
	for(int i=0; i<n; i++) {
		for(int k=0; k<ncov; k++) {
			jnk6 = 0.;
			for (int l=0; l<ncov; l++) jnk6 += omega0inv[k][l]*tempsum2[i][l];
			Wbeta[i*ncov+k] = jnk6;
		}		
	}
	for(int j=0; j<count4; j++) {
		jnk5 = 0.;
		for (int i=0; i<n; i++) jnk5 += pow(W_lambda[i][j],2); 
		lambda10sd[j] = sqrt(jnk5);
	}	
	
	//memory cleanup
	for(int i=0; i<n; i++) delete[] Wc[i];
	delete[] Wc;
	for(int i=0; i<n; i++) delete[] tmp_2[i];
	delete[] tmp_2;	
	for(int i=0; i<n; i++) delete[] W_lambda[i];
	delete[] W_lambda;		
	for(int i=0; i<n; i++) delete[] W_lambda1[i];
	delete[] W_lambda1;	
	for(int i=0; i<n; i++) delete[] W_lambda2[i];
	delete[] W_lambda2;	
	for(int i=0; i<n; i++) delete[] wdm1[i];
	delete[] wdm1;
	for(int i=0; i<n; i++) delete[] yci[i];
	delete[] yci;
	for(int i=0; i<n; i++) delete[] dmci[i];
	delete[] dmci;
	for(int i=0; i<ncov; i++) delete[] omegasum[i];
	delete[] omegasum;
	
	for(int i=0; i<n; i++) {
		for(int j=0; j<njp; j++) delete[] h[i][j];
		delete[] h[i];
	}
	delete[] h;
	for(int i=0; i<n; i++) {
		for(int j=0; j<njp; j++) delete[] jnk0[i][j];
		delete[] jnk0[i];
	}
	delete[] jnk0;
}






