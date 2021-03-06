// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "RcppArmadillo.h"
#include <RcppParallel.h>
#include <Rcpp.h>
#include <math.h>
using namespace RcppParallel;
using namespace Rcpp;
using namespace arma;
using namespace std;
// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]

#include <cmath>
#include <algorithm>



// helper function for taking the average of two numbers
inline double average(double val1, double val2) {
   return (val1 + val2) / 2;
}

  

// helper function for Gibbs sampling

inline double GibbsMCMC(RVector<double> nn, RMatrix<double> data, RVector<double> nnp, RMatrix<double> priordata, RVector<double> priorweight, RMatrix<double> thetaboot,
	RVector<double> bootmean0, RVector<double> bootmean1, RMatrix<double> databoot, RVector<double> scheduleLen, RMatrix<double> priorSched,
	RVector<double> alpha, RVector<double> M_samp, RVector<double> w, std::size_t i) {
   	

	int M = int(M_samp[0]);
	int n = int(nn[0]);
	int np = int(nnp[0]);
	int sN = int(scheduleLen[0]);
   	NumericVector prop0(1,0.0);
   	NumericVector prop1(1,0.0);
	double cov_ind = 0.0;
	NumericVector YIboot(1,0.0);
   	NumericVector datamin(1,0.0);
	NumericVector datamax(1,0.0);
   	NumericVector theta0old(1,0.0);
	NumericVector theta0new(1,0.0);
	NumericVector theta1old(1,0.0);
	NumericVector theta1new(1,0.0);
	NumericVector loglikdiff(1,0.0);
	NumericVector r(1,0.0);
	NumericVector uu(1,0.0);
	NumericVector vv(1,0.0);
	NumericVector postsamples0(M,0.0);
	NumericVector postsamples1(M,0.0);
	NumericVector logpost(M,0.0);
	NumericVector l0(1,0.0);
	NumericVector l1(1,0.0);
	NumericVector u0(1,0.0);
	NumericVector u1(1,0.0);
	NumericVector n0(1,0.0);
	NumericVector n1(1,0.0);
	NumericVector n2(1,0.0);
	NumericVector F0_c0old(1,0.0);
	NumericVector F0_c0new(1,0.0);
	NumericVector F1_c0old(1,0.0);
	NumericVector F1_c0new(1,0.0);
	NumericVector F1_c1old(1,0.0);
	NumericVector F1_c1new(1,0.0);
	NumericVector F2_c1old(1,0.0);
	NumericVector F2_c1new(1,0.0);
	NumericVector n0p(1,0.0);
	NumericVector n1p(1,0.0);
	NumericVector n2p(1,0.0);
	NumericVector F0_c0oldp(1,0.0);
	NumericVector F0_c0newp(1,0.0);
	NumericVector F1_c0oldp(1,0.0);
	NumericVector F1_c0newp(1,0.0);
	NumericVector F1_c1oldp(1,0.0);
	NumericVector F1_c1newp(1,0.0);
	NumericVector F2_c1oldp(1,0.0);
	NumericVector F2_c1newp(1,0.0);
	theta0old(0) = thetaboot(i,0);
	theta1old(0) = thetaboot(i,1);
	NumericVector YI(M,0.0);
	NumericVector YIl(1,0.0);
	NumericVector YIu(1,0.0);
	datamin(0) = databoot(0,2*i+1);
	datamax(0) = databoot(0,2*i+1);
	
	for(int j=0; j<sN; j++){
		if(w[0]<=priorSched(j,0)){
			prop0(0) = priorSched(j,1);	
		}
		if(w[0]<=priorSched(j,2)){
			prop1(0) = priorSched(j,3);	
		}		
	}

	
	for(int k=0; k<n; k++){
		if(databoot(k,2*i)==1){
			n0(0) = n0(0) + 1;
			if(databoot(k,2*i+1)<=theta0old(0)){
				F0_c0old(0) = 	F0_c0old(0)+1;
			}	
		}
		else if(databoot(k,2*i)==2){
			n1(0) = n1(0) + 1;
			if(databoot(k,2*i+1)<=theta0old(0)){
				F1_c0old(0) = 	F1_c0old(0)+1;
			}
			if(databoot(k,2*i+1)<=theta1old(0)){
				F1_c1old(0) = 	F1_c1old(0)+1;
			}	
		}
		else {
			n2(0) = n2(0) + 1;
			if(databoot(k,2*i+1)<=theta1old(0)){
				F2_c1old(0) = 	F2_c1old(0)+1;
			}
		}
		if(datamin(0) > databoot(k,2*i+1)){
			datamin(0) = databoot(k,2*i+1);	
		}
		if(datamax(0) < databoot(k,2*i+1)){
			datamax(0) = databoot(k,2*i+1);	
		}
		YIboot(0) = YIboot(0)+thetaboot(k,2);
	}
	YIboot(0) = YIboot(0)/n;
	F0_c0old(0) = 	F0_c0old(0)/n0(0);
	F1_c0old(0) = 	F1_c0old(0)/n1(0);
	F1_c1old(0) = 	F1_c1old(0)/n1(0);
	F2_c1old(0) = 	F2_c1old(0)/n2(0);
	
	for(int k=0; k<np; k++){
		if(priordata(k,0)==1){
			n0p(0) = n0p(0) + 1;
			if(priordata(k,1)<=theta0old(0)){
				F0_c0oldp(0) = 	F0_c0oldp(0)+1;
			}	
		}
		else if(priordata(k,0)==2){
			n1p(0) = n1p(0) + 1;
			if(priordata(k,1)<=theta0old(0)){
				F1_c0oldp(0) = 	F1_c0oldp(0)+1;
			}
			if(priordata(k,1)<=theta1old(0)){
				F1_c1oldp(0) = 	F1_c1oldp(0)+1;
			}	
		}
		else {
			n2p(0) = n2p(0) + 1;
			if(priordata(k,1)<=theta1old(0)){
				F2_c1oldp(0) = 	F2_c1oldp(0)+1;
			}
		}
	}
	F0_c0oldp(0) = 	F0_c0oldp(0)/n0p(0);
	F1_c0oldp(0) = 	F1_c0oldp(0)/n1p(0);
	F1_c1oldp(0) = 	F1_c1oldp(0)/n1p(0);
	F2_c1oldp(0) = 	F2_c1oldp(0)/n2p(0);

	
	for(int j=0; j<M; j++) {
		theta0new(0) = R::rnorm(theta0old(0), prop0(0));
		theta1new(0) = R::rnorm(theta1old(0), prop1(0));
	
		if( (theta0new(0)>theta1new(0)) ||  (theta0new(0)<datamin(0)) || (theta1new(0)>datamax(0))   ){
			postsamples0(j) = theta0old(0);
			postsamples1(j) = theta1old(0);
			logpost(j) = w[0]*n*(F0_c0old(0)+F1_c1old(0)-F1_c0old(0)-F2_c1old(0))+w[0]*priorweight[0]*np*(F0_c0oldp(0)+F1_c1oldp(0)-F1_c0oldp(0)-F2_c1oldp(0));
			YI(j) = F0_c0old(0) + F1_c1old(0) - F1_c0old(0) - F2_c1old(0);
		}else {
		
		loglikdiff(0) = 0.0;
		for(int k=0; k<n; k++){
			if(databoot(k,2*i)==1){
				if(databoot(k,2*i+1)<=theta0new(0)){
					F0_c0new(0) = 	F0_c0new(0)+1.0;
				}	
			}
			else if(databoot(k,2*i)==2){
				if(databoot(k,2*i+1)<=theta0new(0)){
					F1_c0new(0) = 	F1_c0new(0)+1.0;
				}
				if(databoot(k,2*i+1)<=theta1new(0)){
					F1_c1new(0) = 	F1_c1new(0)+1.0;
				}	
			}
			else {
				if(databoot(k,2*i+1)<=theta1new(0)){
					F2_c1new(0) = 	F2_c1new(0)+1.0;
				}
			}
		}
		F0_c0new(0) = 	F0_c0new(0)/n0(0);
		F1_c0new(0) = 	F1_c0new(0)/n1(0);
		F1_c1new(0) = 	F1_c1new(0)/n1(0);
		F2_c1new(0) = 	F2_c1new(0)/n2(0);
		
		for(int k=0; k<np; k++){
			if(priordata(k,0)==1){
				if(priordata(k,1)<=theta0new(0)){
					F0_c0newp(0) = 	F0_c0newp(0)+1.0;
				}	
			}
			else if(priordata(k,0)==2){
				if(priordata(k,1)<=theta0new(0)){
					F1_c0newp(0) = 	F1_c0newp(0)+1.0;
				}
				if(priordata(k,1)<=theta1new(0)){
					F1_c1newp(0) = 	F1_c1newp(0)+1.0;
				}	
			}
			else {
				if(priordata(k,1)<=theta1new(0)){
					F2_c1newp(0) = 	F2_c1newp(0)+1.0;
				}
			}
		}
		F0_c0newp(0) = 	F0_c0newp(0)/n0p(0);
		F1_c0newp(0) = 	F1_c0newp(0)/n1p(0);
		F1_c1newp(0) = 	F1_c1newp(0)/n1p(0);
		F2_c1newp(0) = 	F2_c1newp(0)/n2p(0);
	
		loglikdiff(0) = w[0]*n*((F0_c0new(0)+F1_c1new(0)-F1_c0new(0)-F2_c1new(0))-(F0_c0old(0)+F1_c1old(0)-F1_c0old(0)-F2_c1old(0)))  +   w[0]*priorweight[0]*np*((F0_c0newp(0)+F1_c1newp(0)-F1_c0newp(0)-F2_c1newp(0))-(F0_c0oldp(0)+F1_c1oldp(0)-F1_c0oldp(0)-F2_c1oldp(0)));
		loglikdiff(0) = fmin(std::exp(loglikdiff(0)), 1.0);
		uu[0] = R::runif(0.0,1.0);
		if(uu(0) <= loglikdiff(0)) {
			postsamples0(j) = theta0new(0);
			postsamples1(j) = theta1new(0);
			logpost(j) = w[0]*n*(F0_c0new(0)+F1_c1new(0)-F1_c0new(0)-F2_c1new(0))+w[0]*priorweight[0]*np*(F0_c0newp(0)+F1_c1newp(0)-F1_c0newp(0)-F2_c1newp(0));
			theta0old(0) = theta0new(0);
			theta1old(0) = theta1new(0);
			F0_c0old(0)=F0_c0new(0);
			F1_c0old(0)=F1_c0new(0);
			F1_c1old(0)=F1_c1new(0);
			F2_c1old(0)=F2_c1new(0);
			F0_c0oldp(0)=F0_c0newp(0);
			F1_c0oldp(0)=F1_c0newp(0);
			F1_c1oldp(0)=F1_c1newp(0);
			F2_c1oldp(0)=F2_c1newp(0);
			YI(j) = F0_c0new(0) + F1_c1new(0) - F1_c0new(0) - F2_c1new(0);
		}
		else {
			postsamples0(j) = theta0old(0);
			postsamples1(j) = theta1old(0);
			logpost(j) = w[0]*n*(F0_c0old(0)+F1_c1old(0)-F1_c0old(0)-F2_c1old(0))+w[0]*priorweight[0]*np*(F0_c0oldp(0)+F1_c1oldp(0)-F1_c0oldp(0)-F2_c1oldp(0));
			YI(j) = F0_c0old(0) + F1_c1old(0) - F1_c0old(0) - F2_c1old(0);
		}
		}
	}	
/*
	double templogpost;
	double tempYI;
	double temppost0;
	double temppost1;
	bool swapped;	
        for (int j = 0; j < M-1; j++){ 
		 swapped = false; 
       		 for (int k = 0; k < M-j-1; k++){  
           		 if (logpost(k) > logpost(k+1)){ 
				 swapped = true;
				 templogpost = logpost(k);
				 logpost(k) = logpost(k+1);
				 logpost(k+1) = templogpost;
				 tempYI = YI(k);
				 YI(k) = YI(k+1);
				 YI(k+1) = tempYI;
				 temppost0 = postsamples0(k);
				 postsamples0(k) = postsamples0(k+1);
				 postsamples0(k+1) = temppost0;				 
				 temppost1 = postsamples1(k);
				 postsamples1(k) = postsamples1(k+1);
				 postsamples1(k+1) = temppost1;
			 }
		 }
		 if (swapped == false) 
			break;
        }
			
	l0[0] = 100000;
	u0[0] = -100000;
	l1[0] = 100000;
	u1[0] = -100000;
	YIl[0] = 100000;
	YIu[0] = -100000;
	for(int i = (0.05*M-1); i<M; i++){
		if(l0[0]>postsamples0(i)){
			l0[0]=postsamples0(i);
		}
		if(u0[0]<postsamples0(i)){
			u0[0]=postsamples0(i);
		}
		if(l1[0]>postsamples1(i)){
			l1[0]=postsamples1(i);
		}
		if(u1[0]<postsamples1(i)){
			u1[0]=postsamples1(i);
		}
		if(YIl[0]>YI(i)){
			YIl[0]=YI(i);
		}
		if(YIu[0]<YI(i)){
			YIu[0]=YI(i);
		}
	}*/
	
	std::sort(postsamples0.begin(), postsamples0.end());
	std::sort(postsamples1.begin(), postsamples1.end());
	std::sort(YI.begin(), YI.end());
	l0[0] = postsamples0(M*.025-1);
	u0[0] = postsamples0(M*.975-1);
	l1[0] = postsamples1(M*.025-1);
	u1[0] = postsamples1(M*.975-1);
	YIl[0] = YI[M*.025-1];
	YIu[0] = YI[M*.975-1];
	/*if ( (YIl[0] < YIboot[0]) && (YIu[0] > YIboot[0]) ){
			cov_ind = 1.0;
	} else {cov_ind = 0.0;}*/
	if ( (l0[0] < bootmean0[0]) && (u0[0] > bootmean0[0]) ){
			cov_ind = 1.0;
	} else {cov_ind = 0.0;}
	return cov_ind;

	
}


inline double GibbsMCMCkde(RVector<double> nn, RMatrix<double> data, RVector<double> nnp, RMatrix<double> priordata, RVector<double> priorweight, RMatrix<double> thetaboot,
	RVector<double> bootmean0, RVector<double> bootmean1, RVector<double> bootmeanYI, RVector<double> kdecdflen, RMatrix<double> kdecdfboot1, RMatrix<double> kdecdfboot2, RMatrix<double> kdecdfboot3, RMatrix<double> kdecdfboot1p, RMatrix<double> kdecdfboot2p, RMatrix<double> kdecdfboot3p, RVector<double> scheduleLen, RMatrix<double> priorSched,
	RVector<double> alpha, RVector<double> M_samp, RVector<double> w, std::size_t i) {
   	

	int M = int(M_samp[0]);
	int n = int(nn[0]);
	int np = int(nnp[0]);
	int sN = int(scheduleLen[0]);
	int kdeN = int(kdecdflen[0]);
   	NumericVector prop0(1,0.0);
   	NumericVector prop1(1,0.0);
	double cov_ind = 0.0;
   	NumericVector theta0old(1,0.0);
	NumericVector theta0new(1,0.0);
	NumericVector theta1old(1,0.0);
	NumericVector theta1new(1,0.0);
	NumericVector loglikdiff(1,0.0);
	NumericVector r(1,0.0);
	NumericVector uu(1,0.0);
	NumericVector vv(1,0.0);
	NumericVector postsamples0(M,0.0);
	NumericVector postsamples1(M,0.0);
	NumericVector logpost(M,0.0);
	NumericVector l0(1,0.0);
	NumericVector l1(1,0.0);
	NumericVector u0(1,0.0);
	NumericVector u1(1,0.0);
	NumericVector F0_c0old(1,0.0);
	NumericVector F0_c0new(1,0.0);
	NumericVector F1_c0old(1,0.0);
	NumericVector F1_c0new(1,0.0);
	NumericVector F1_c1old(1,0.0);
	NumericVector F1_c1new(1,0.0);
	NumericVector F2_c1old(1,0.0);
	NumericVector F2_c1new(1,0.0);
	NumericVector F0_c0oldp(1,0.0);
	NumericVector F0_c0newp(1,0.0);
	NumericVector F1_c0oldp(1,0.0);
	NumericVector F1_c0newp(1,0.0);
	NumericVector F1_c1oldp(1,0.0);
	NumericVector F1_c1newp(1,0.0);
	NumericVector F2_c1oldp(1,0.0);
	NumericVector F2_c1newp(1,0.0);
	theta0old(0) = thetaboot(i,0);
	theta1old(0) = thetaboot(i,1);
	NumericVector YI(M,0.0);
	NumericVector YIl(1,0.0);
	NumericVector YIu(1,0.0);

	
	for(int j=0; j<sN; j++){
		if(w[0]<=priorSched(j,0)){
			prop0(0) = priorSched(j,1);	
		}
		if(w[0]<=priorSched(j,2)){
			prop1(0) = priorSched(j,3);	
		}		
	}

	
	int k = 0;
	while(k < kdeN){
		if(theta0old(0)<=kdecdfboot1(k,2*i)){
			F0_c0old(0) = kdecdfboot1(k,2*i+1);
			break;
		}	
		k = k + 1;
	}
	k = 0;
	while(k < kdeN){
		if(theta0old(0)<=kdecdfboot2(k,2*i)){
			F1_c0old(0) = kdecdfboot2(k,2*i+1);
			break;
		}
		k = k + 1;
	}
	k = 0;
	while(k < kdeN){
		if(theta1old(0)<=kdecdfboot2(k,2*i)){
			F1_c1old(0) = kdecdfboot2(k,2*i+1);
			break;
		}
		k = k + 1;
	}
	k = 0;
	while(k < kdeN){
		if(theta1old(0)<=kdecdfboot3(k,2*i)){
			F2_c1old(0) = kdecdfboot3(k,2*i+1);
			break;
		}	
		k = k + 1;
	}	
	
	
	if(priorweight[0]>0.0){
		k = 0;
		while(k < kdeN){
			if(theta0old(0)<=kdecdfboot1p(k,2*i)){
				F0_c0oldp(0) = kdecdfboot1p(k,2*i+1);
				break;
			}
			k = k + 1;
		}
		k = 0;
		while(k < kdeN){
			if(theta0old(0)<=kdecdfboot2p(k,2*i)){
				F1_c0oldp(0) = kdecdfboot2p(k,2*i+1);
				break;
			}
			k = k + 1;
		}
		k = 0;
		while(k < kdeN){
			if(theta1old(0)<=kdecdfboot2p(k,2*i)){
				F1_c1oldp(0) = kdecdfboot2p(k,2*i+1);
				break;
			}
			k = k + 1;
		}
		k = 0;
		while(k < kdeN){
			if(theta1old(0)<=kdecdfboot3p(k,2*i)){
				F2_c1oldp(0) = kdecdfboot3p(k,2*i+1);
				break;
			}
			k = k + 1;
		}	
	}


	
	for(int j=0; j<M; j++) {
		theta0new(0) = R::rnorm(theta0old(0), prop0(0));
		theta1new(0) = R::rnorm(theta1old(0), prop1(0));
	
		if( (theta0new(0)>theta1new(0)) || (theta0new(0) < kdecdfboot1(0,2*i)) || (theta0new(0) > kdecdfboot3(kdeN-1,2*i))){
			postsamples0(j) = theta0old(0);
			postsamples1(j) = theta1old(0);
			logpost(j) = w[0]*n*(F0_c0old(0)+F1_c1old(0)-F1_c0old(0)-F2_c1old(0))+w[0]*priorweight[0]*np*(F0_c0oldp(0)+F1_c1oldp(0)-F1_c0oldp(0)-F2_c1oldp(0));
			YI(j) = F0_c0old(0) + F1_c1old(0) - F1_c0old(0) - F2_c1old(0);
		}else {
		
		loglikdiff(0) = 0.0;
		k = 0;
		while(k < kdeN){
			if(theta0new(0)<=kdecdfboot1(k,2*i)){
				F0_c0new(0) = kdecdfboot1(k,2*i+1);
				break;
			}
			k = k + 1;
		}
		k = 0;
		while(k < kdeN){
			if(theta0new(0)<=kdecdfboot2(k,2*i)){
				F1_c0new(0) = kdecdfboot2(k,2*i+1);
				break;
			}	
			k = k + 1;
		}
		k = 0;
		while(k < kdeN){
			if(theta1new(0)<=kdecdfboot2(k,2*i)){
				F1_c1new(0) = kdecdfboot2(k,2*i+1);
				break;
			}
			k = k + 1;
		}
		k = 0;
		while(k < kdeN){
			if(theta1new(0)<=kdecdfboot3(k,2*i)){
				F2_c1new(0) = kdecdfboot3(k,2*i+1);
				break;
			}
			k = k + 1;
		}	
	
	
		if(priorweight[0]>0.0){
			k = 0;
			while(k < kdeN){
				if(theta0new(0)<=kdecdfboot1p(k,2*i)){
					F0_c0newp(0) = kdecdfboot1p(k,2*i+1);
					break;
				}
				k = k + 1;
			}
			k = 0;
			while(k < kdeN){
				if(theta0new(0)<=kdecdfboot2p(k,2*i)){
					F1_c0newp(0) = kdecdfboot2p(k,2*i+1);
					break;
				}
				k = k + 1;
			}
			k = 0;
			while(k < kdeN){
				if(theta1new(0)<=kdecdfboot2p(k,2*i)){
					F1_c1newp(0) = kdecdfboot2p(k,2*i+1);
					break;
				}
				k = k + 1;
			}
			k = 0;
			while(k < kdeN){
				if(theta1new(0)<=kdecdfboot3p(k,2*i)){
					F2_c1newp(0) = kdecdfboot3p(k,2*i+1);
					break;
				}
				k = k + 1;
			}	
		}
	
		loglikdiff(0) = w[0]*n*((F0_c0new(0)+F1_c1new(0)-F1_c0new(0)-F2_c1new(0))-(F0_c0old(0)+F1_c1old(0)-F1_c0old(0)-F2_c1old(0)))  +   w[0]*priorweight[0]*np*((F0_c0newp(0)+F1_c1newp(0)-F1_c0newp(0)-F2_c1newp(0))-(F0_c0oldp(0)+F1_c1oldp(0)-F1_c0oldp(0)-F2_c1oldp(0)));
		loglikdiff(0) = fmin(std::exp(loglikdiff(0)), 1.0);
		uu[0] = R::runif(0.0,1.0);
		if(uu(0) <= loglikdiff(0)) {
			postsamples0(j) = theta0new(0);
			postsamples1(j) = theta1new(0);
			logpost(j) = w[0]*n*(F0_c0new(0)+F1_c1new(0)-F1_c0new(0)-F2_c1new(0))+w[0]*priorweight[0]*np*(F0_c0newp(0)+F1_c1newp(0)-F1_c0newp(0)-F2_c1newp(0));
			theta0old(0) = theta0new(0);
			theta1old(0) = theta1new(0);
			F0_c0old(0)=F0_c0new(0);
			F1_c0old(0)=F1_c0new(0);
			F1_c1old(0)=F1_c1new(0);
			F2_c1old(0)=F2_c1new(0);
			F0_c0oldp(0)=F0_c0newp(0);
			F1_c0oldp(0)=F1_c0newp(0);
			F1_c1oldp(0)=F1_c1newp(0);
			F2_c1oldp(0)=F2_c1newp(0);
			YI(j) = F0_c0new(0) + F1_c1new(0) - F1_c0new(0) - F2_c1new(0);
		}
		else {
			postsamples0(j) = theta0old(0);
			postsamples1(j) = theta1old(0);
			logpost(j) = w[0]*n*(F0_c0old(0)+F1_c1old(0)-F1_c0old(0)-F2_c1old(0))+w[0]*priorweight[0]*np*(F0_c0oldp(0)+F1_c1oldp(0)-F1_c0oldp(0)-F2_c1oldp(0));
			YI(j) = F0_c0old(0) + F1_c1old(0) - F1_c0old(0) - F2_c1old(0);
		}
		}
	}	
/*
	double templogpost;
	double tempYI;
	double temppost0;
	double temppost1;
	bool swapped;	
        for (int j = 0; j < M-1; j++){ 
		 swapped = false; 
       		 for (int k = 0; k < M-j-1; k++){  
           		 if (logpost(k) > logpost(k+1)){ 
				 swapped = true;
				 templogpost = logpost(k);
				 logpost(k) = logpost(k+1);
				 logpost(k+1) = templogpost;
				 tempYI = YI(k);
				 YI(k) = YI(k+1);
				 YI(k+1) = tempYI;
				 temppost0 = postsamples0(k);
				 postsamples0(k) = postsamples0(k+1);
				 postsamples0(k+1) = temppost0;				 
				 temppost1 = postsamples1(k);
				 postsamples1(k) = postsamples1(k+1);
				 postsamples1(k+1) = temppost1;
			 }
		 }
		 if (swapped == false) 
			break;
        }
			
	l0[0] = 100000;
	u0[0] = -100000;
	l1[0] = 100000;
	u1[0] = -100000;
	YIl[0] = 100000;
	YIu[0] = -100000;
	for(int i = (0.05*M-1); i<M; i++){
		if(l0[0]>postsamples0(i)){
			l0[0]=postsamples0(i);
		}
		if(u0[0]<postsamples0(i)){
			u0[0]=postsamples0(i);
		}
		if(l1[0]>postsamples1(i)){
			l1[0]=postsamples1(i);
		}
		if(u1[0]<postsamples1(i)){
			u1[0]=postsamples1(i);
		}
		if(YIl[0]>YI(i)){
			YIl[0]=YI(i);
		}
		if(YIu[0]<YI(i)){
			YIu[0]=YI(i);
		}
	}*/
	
	std::sort(postsamples0.begin(), postsamples0.end());
	std::sort(postsamples1.begin(), postsamples1.end());
	std::sort(YI.begin(), YI.end());
	l0[0] = postsamples0(M*.025-1);
	u0[0] = postsamples0(M*.975-1);
	l1[0] = postsamples1(M*.025-1);
	u1[0] = postsamples1(M*.975-1);
	YIl[0] = YI[M*.025-1];
	YIu[0] = YI[M*.975-1];
	/*if ( (YIl[0] < bootmeanYI[0]) && (YIu[0] > bootmeanYI[0]) ){
			cov_ind = 1.0;
	} else {cov_ind = 0.0;}*/
	if ( (l0[0] < bootmean0[0]) && (u0[0] > bootmean0[0]) ){
			cov_ind = 1.0;
	} else {cov_ind = 0.0;}
	return cov_ind;

	
}



// Proposal Schedule as input
// [[Rcpp::export]]
Rcpp::List GridSearchKDE(int n, NumericVector mesh1, NumericVector mesh2, NumericVector mesh3, NumericVector cdf1, NumericVector cdf2, NumericVector cdf3, int n12, int n23, NumericVector d12, NumericVector d23 ) {
   	
	List result;
   	NumericVector YIhat(1,0.0);	
	NumericVector theta0hat(1,0.0);
	NumericVector theta1hat(1,0.0);
	NumericVector YI(1,0.0);	
	NumericVector theta0(1,0.0);
	NumericVector theta1(1,0.0);
	NumericVector diff1(1,10.0);
	NumericVector diff2a(1,10.0);
	NumericVector diff2b(1,10.0);
	NumericVector diff3(1,10.0);
	NumericVector diff(1,0.0);
	int index1 = 0;
	int index2a = 0;
	int index2b = 0;
	int index3 = 0;
	
	for(int i = 0; i < n12; i++){
		for(int j = 0; j < n23; j++){
			theta0(0) = d12(i);
			theta1(0) = d23(j);
			if(theta0(0)<theta1(0)){
				for(int k = 0; k < n; k++){
					diff(0) = abs(mesh1(k) - theta0(0)); 
					if(diff(0)<diff1(0)){
						index1 = k;	
						diff1(0) = diff(0);
					}
					diff(0) = abs(mesh2(k) - theta0(0)); 
					if(diff(0)<diff2a(0)){
						index2a = k;	
						diff2a(0) = diff(0);
					}
					diff(0) = abs(mesh2(k) - theta1(0)); 
					if(diff(0)<diff2b(0)){
						index2b = k;	
						diff2b(0) = diff(0);
					}
					diff(0) = abs(mesh3(k) - theta1(0)); 
					if(diff(0)<diff3(0)){
						index3 = k;	
						diff3(0) = diff(0);
					}
				}
				YI(0) = cdf1(index1)-cdf2(index2a)+cdf2(index2b)-cdf3(index3);
				if(YI(0)>YIhat(0)){
					YIhat(0) = YI(0);
					theta0hat(0) = theta0(0);
					theta1hat(0) = theta1(0);
				}	
			}
		}
	}
	
	result = Rcpp::List::create(Rcpp::Named("YI") = YIhat, Rcpp::Named("theta0") = theta0hat, Rcpp::Named("theta1") = theta1hat);

	return result;
	
}
   	



inline std::vector<double> GibbsMCMCsmooth(RVector<double> nn, RMatrix<double> data1, RMatrix<double> data2, RMatrix<double> thetaboot,
	RVector<double> bootmean0, RVector<double> bootmean1, RMatrix<double> databoot1, RMatrix<double> databoot2, RVector<double> normprior, 
	RVector<double> scheduleLen, RMatrix<double> propSched, RVector<double> ddelta, RVector<double> alpha, RVector<double> M_samp, 
	RVector<double> w, std::size_t i) {
   	
	double pi = 3.14159265358979323846;
	int M = int(M_samp[0]);
	int n1 = int(nn[0]);
	int n2 = int(nn[1]);
	int sN = int(scheduleLen[0]);
	double delta0 = double(ddelta[0]);
	double delta1 = double(ddelta[1]);
	double w0 = double(w[0]);
	double w1 = double(w[1]);
   	NumericVector prop0(1,0.0);
   	NumericVector prop1(1,0.0);
	std::vector<double> cov_ind;
	double z1 = 0.0; double z2 = 0.0;
   	NumericVector theta0old(1,0.0);
	NumericVector theta0new(1,0.0);
	NumericVector theta1old(1,0.0);
	NumericVector theta1new(1,0.0);
	NumericVector theta0temp(1,0.0);
	NumericVector theta1temp(1,0.0);
	NumericVector loglikdiff(1,0.0);
	NumericVector loss1temp(1,0.0);
	NumericVector loss2temp(1,0.0);	
	NumericVector loss1old(1,0.0);
	NumericVector loss1new(1,0.0);
	NumericVector loss2old(1,0.0);
	NumericVector loss2new(1,0.0);
	NumericVector r(1,0.0);
	NumericVector uu(1,0.0);
	NumericVector postsamples0(M,0.0);
	NumericVector postsamples1(M,0.0);
	NumericVector l0(1,0.0);
	NumericVector l1(1,0.0);
	NumericVector u0(1,0.0);
	NumericVector u1(1,0.0);
	theta0old(0) = thetaboot(i,0);
	theta1old(0) = thetaboot(i,1);



	
	for(int j=0; j<sN; j++){
		if(w0<=propSched(j,0)){
			prop0(0) = propSched(j,1);	
		}
		if(w1<=propSched(j,2)){
			prop1(0) = propSched(j,3);	
		}		
	}
	if(w0>1.0){
		prop0(0) = propSched(0,1);	
	}
	if(w1>1.0){
		prop1(0) = propSched(0,3);	
	}


	for(int k=0; k<n1; k++){
		z1 = databoot1(k,2*i)*(databoot1(k,2*i+1)-theta0old(0));
		loss1temp(0)=0;
		if((z1>-delta0/2) && (z1<=delta0/2)){
			loss1temp(0) = 0.5 + 0.5*sin((pi/delta0)*(z1)+pi);
		}else if(z1<0){
			loss1temp(0) = 1.0;	
		}
		loss1old(0) = loss1old(0) + loss1temp(0);
	}
	for(int k=0; k<n2; k++){		
		z2 = databoot2(k,2*i)*(databoot2(k,2*i+1)-theta1old(0));
		loss2temp(0)=0;
		if((z2>-delta1/2) && (z2<=delta1/2)){
			loss2temp(0) = 0.5 + 0.5*sin((pi/delta1)*(z2)+pi);
		}else if(z2<0){
			loss2temp(0) = 1.0;	
		}
		loss2old(0) = loss2old(0) + loss2temp(0);
	}

	for(int j=0; j<M; j++) {
		theta0temp(0) = R::rnorm(theta0old(0), prop0(0));
		theta1temp(0) = R::rnorm(theta1old(0), prop1(0));
		theta0new(0) = min(theta0temp(0), theta1temp(0));
		theta1new(0) = max(theta0temp(0), theta1temp(0));
		
		loss1new(0) = 0.0;
		loss2new(0) = 0.0;
		loglikdiff(0) = 0.0;
		for(int k=0; k<n1; k++){
			z1 = databoot1(k,2*i)*(databoot1(k,2*i+1)-theta0new(0));
			loss1temp(0)=0;
			if((z1>-delta0/2) && (z1<=delta0/2)){
				loss1temp(0) = 0.5 + 0.5*sin((pi/delta0)*(z1)+pi);
			}else if(z1<0){
				loss1temp(0) = 1.0;	
			}
			loss1new(0) = loss1new(0) + loss1temp(0);
		}
		for(int k=0; k<n2; k++){		
			z2 = databoot2(k,2*i)*(databoot2(k,2*i+1)-theta1new(0));
			loss2temp(0)=0;
			if((z2>-delta1/2) && (z2<=delta1/2)){
				loss2temp(0) = 0.5 + 0.5*sin((pi/delta1)*(z2)+pi);
			}else if(z2<0){
				loss2temp(0) = 1.0;	
			}
			loss2new(0) = loss2new(0) + loss2temp(0);
		}

		loglikdiff(0) = -w1*(loss2new(0)-loss2old(0))-w0*(loss1new(0)-loss1old(0));
		loglikdiff(0) = fmin(std::exp(loglikdiff(0))*((R::dnorm(theta0new(0),normprior[0],normprior[1],0)*R::dnorm(theta1new(0),normprior[2],normprior[3],0)+R::dnorm(theta1new(0),normprior[0],normprior[1],0)*R::dnorm(theta0new(0),normprior[2],normprior[3],0))/(R::dnorm(theta0old(0),normprior[0],normprior[1],0)*R::dnorm(theta1old(0),normprior[2],normprior[3],0)+R::dnorm(theta1old(0),normprior[0],normprior[1],0)*R::dnorm(theta0old(0),normprior[2],normprior[3],0))), 1.0);
		uu[0] = R::runif(0.0,1.0);
		if(uu(0) <= loglikdiff(0)) {
			postsamples0(j) = theta0new(0);
			postsamples1(j) = theta1new(0);
			theta0old(0) = theta0new(0);
			theta1old(0) = theta1new(0);
			loss2old(0) = loss2new(0);
			loss1old(0) = loss1new(0);	
		}
		else {
			postsamples0(j) = theta0old(0);
			postsamples1(j) = theta1old(0);
		}
	}	
/*
	double templogpost;
	double tempYI;
	double temppost0;
	double temppost1;
	bool swapped;	
        for (int j = 0; j < M-1; j++){ 
		 swapped = false; 
       		 for (int k = 0; k < M-j-1; k++){  
           		 if (logpost(k) > logpost(k+1)){ 
				 swapped = true;
				 templogpost = logpost(k);
				 logpost(k) = logpost(k+1);
				 logpost(k+1) = templogpost;
				 tempYI = YI(k);
				 YI(k) = YI(k+1);
				 YI(k+1) = tempYI;
				 temppost0 = postsamples0(k);
				 postsamples0(k) = postsamples0(k+1);
				 postsamples0(k+1) = temppost0;				 
				 temppost1 = postsamples1(k);
				 postsamples1(k) = postsamples1(k+1);
				 postsamples1(k+1) = temppost1;
			 }
		 }
		 if (swapped == false) 
			break;
        }
			
	l0[0] = 100000;
	u0[0] = -100000;
	l1[0] = 100000;
	u1[0] = -100000;
	YIl[0] = 100000;
	YIu[0] = -100000;
	for(int i = (0.05*M-1); i<M; i++){
		if(l0[0]>postsamples0(i)){
			l0[0]=postsamples0(i);
		}
		if(u0[0]<postsamples0(i)){
			u0[0]=postsamples0(i);
		}
		if(l1[0]>postsamples1(i)){
			l1[0]=postsamples1(i);
		}
		if(u1[0]<postsamples1(i)){
			u1[0]=postsamples1(i);
		}
		if(YIl[0]>YI(i)){
			YIl[0]=YI(i);
		}
		if(YIu[0]<YI(i)){
			YIu[0]=YI(i);
		}
	}*/
	
	std::sort(postsamples0.begin(), postsamples0.end());
	std::sort(postsamples1.begin(), postsamples1.end());
	l0[0] = postsamples0(M*.025-1);
	u0[0] = postsamples0(M*.975-1);
	l1[0] = postsamples1(M*.025-1);
	u1[0] = postsamples1(M*.975-1);
	if ( (l0[0] < bootmean0[0]) && (u0[0] > bootmean0[0]) ){
		cov_ind.push_back(1.0);
	} else {cov_ind.push_back(0.0);}
	if ( (l1[0] < bootmean1[0]) && (u1[0] > bootmean1[0]) ){
		cov_ind.push_back(1.0);
	} else {cov_ind.push_back(0.0);}
	
	return cov_ind;	
}




// Proposal Schedule as input
// [[Rcpp::export]]
Rcpp::List GibbsMCMC2smooth(NumericVector nn, NumericMatrix data1, NumericMatrix data2, NumericVector bootmean0, NumericVector bootmean1, NumericVector normprior, NumericVector scheduleLen, NumericMatrix propSched, NumericVector alpha, NumericVector ddelta, NumericVector M_samp, NumericVector w) {
	
	List result;
	double pi = 3.14159265358979323846;
	int M = int(M_samp[0]);
	int n1 = int(nn[0]);
	int n2 = int(nn[1]);
	int sN = int(scheduleLen[0]);
	double delta0 = double(ddelta[0]);
	double delta1 = double(ddelta[1]);
	double w0 = double(w[0]);
	double w1 = double(w[1]);
   	NumericVector prop0(1,0.0);
   	NumericVector prop1(1,0.0);
	double z1 = 0.0; double z2 = 0.0;
   	NumericVector theta0old(1,0.0);
	NumericVector theta0new(1,0.0);
	NumericVector theta1old(1,0.0);
	NumericVector theta1new(1,0.0);
	NumericVector theta0temp(1,0.0);
	NumericVector theta1temp(1,0.0);
	NumericVector loglikdiff(1,0.0);
	NumericVector loglikdiff1(1,0.0);
	NumericVector loss1temp(1,0.0);
	NumericVector loss2temp(1,0.0);	
	NumericVector loss1old(1,0.0);
	NumericVector loss1new(1,0.0);
	NumericVector loss2old(1,0.0);
	NumericVector loss2new(1,0.0);
	NumericVector r(1,0.0);
	NumericVector uu(1,0.0);
	NumericVector postsamples0(M,0.0);
	NumericVector postsamples1(M,0.0);
	NumericVector logpost(M,0.0);
	NumericVector l0(1,0.0);
	NumericVector l1(1,0.0);
	NumericVector u0(1,0.0);
	NumericVector u1(1,0.0);
	theta0old(0) = bootmean0(0);
	theta1old(0) = bootmean1(0);
	NumericVector acc(1, 0.0);



	
	for(int j=0; j<sN; j++){
		if(w0<=propSched(j,0)){
			prop0(0) = propSched(j,1);	
		}
		if(w1<=propSched(j,2)){
			prop1(0) = propSched(j,3);	
		}		
	}
	if(w0>1.0){
		prop0(0) = propSched(0,1);	
	}
	if(w1>1.0){
		prop1(0) = propSched(0,3);	
	}


	for(int k=0; k<n1; k++){
		z1 = data1(k,0)*(data1(k,1)-theta0old(0));
		loss1temp(0)=0;
		if((z1>-delta0/2) && (z1<=delta0/2)){
			loss1temp(0) = 0.5 + 0.5*sin((pi/delta0)*(z1)+pi);
		}else if(z1<0){
			loss1temp(0) = 1.0;	
		}
		loss1old(0) = loss1old(0) + loss1temp(0);
	}
	for(int k=0; k<n2; k++){
		z2 = data2(k,0)*(data2(k,1)-theta1old(0));
		loss2temp(0)=0;
		if((z2>-delta1/2) && (z2<=delta1/2)){
			loss2temp(0) = 0.5 + 0.5*sin((pi/delta1)*(z2)+pi);
		}else if(z2<0){
			loss2temp(0) = 1.0;	
		}
		loss2old(0) = loss2old(0) + loss2temp(0);
	}
	

	for(int j=0; j<M; j++) {
		theta0temp(0) = R::rnorm(theta0old(0), prop0(0));
		theta1temp(0) = R::rnorm(theta1old(0), prop1(0));
		theta0new(0) = min(theta0temp(0), theta1temp(0));
		theta1new(0) = max(theta0temp(0), theta1temp(0));

		
		loss1new(0) = 0.0;
		loss2new(0) = 0.0;
		for(int k=0; k<n1; k++){
			z1 = data1(k,0)*(data1(k,1)-theta0new(0));
			loss1temp(0)=0;
			if((z1>-delta0/2) && (z1<=delta0/2)){
				loss1temp(0) = 0.5 + 0.5*sin((pi/delta0)*(z1)+pi);
			}else if(z1<0){
				loss1temp(0) = 1.0;	
			}
			loss1new(0) = loss1new(0) + loss1temp(0);
		}
		for(int k=0; k<n2; k++){
			z2 = data2(k,0)*(data2(k,1)-theta1new(0));
			loss2temp(0)=0;
			if((z2>-delta1/2) && (z2<=delta1/2)){
				loss2temp(0) = 0.5 + 0.5*sin((pi/delta1)*(z2)+pi);
			}else if(z2<0){
				loss2temp(0) = 1.0;	
			}
			loss2new(0) = loss2new(0) + loss2temp(0);
		}
		loglikdiff(0) = 0.0;
		loglikdiff1(0) = 0.0;
		loglikdiff(0) = -w1*(loss2new(0)-loss2old(0))-w0*(loss1new(0)-loss1old(0));
		loglikdiff1(0) = fmin(std::exp(loglikdiff(0))*((R::dnorm(theta0new(0),normprior[0],normprior[1],false)*R::dnorm(theta1new(0),normprior[2],normprior[3],false)+R::dnorm(theta1new(0),normprior[0],normprior[1],false)*R::dnorm(theta0new(0),normprior[2],normprior[3],false))/(R::dnorm(theta0old(0),normprior[0],normprior[1],false)*R::dnorm(theta1old(0),normprior[2],normprior[3],false)+R::dnorm(theta1old(0),normprior[0],normprior[1],false)*R::dnorm(theta0old(0),normprior[2],normprior[3],false))), 1.0);
		uu[0] = R::runif(0.0,1.0);
		if(uu(0) <= loglikdiff1(0)) {
			postsamples0(j) = theta0new(0);
			postsamples1(j) = theta1new(0);
			logpost(j) = log(loglikdiff(0));
			theta0old(0) = theta0new(0);
			theta1old(0) = theta1new(0);
			loss2old(0) = loss2new(0);
			loss1old(0) = loss1new(0);
			acc(0) = acc(0)+1;
		}
		else {
			postsamples0(j) = theta0old(0);
			postsamples1(j) = theta1old(0);
			if(j>0){
				logpost(j) = logpost(j-1);
			}else {
				logpost(j) = logpost(0);	
			}
		}
	}	

	std::sort(postsamples0.begin(), postsamples0.end());
	std::sort(postsamples1.begin(), postsamples1.end());
	l0[0] = postsamples0(M*.025-1);
	u0[0] = postsamples0(M*.975-1);
	l1[0] = postsamples1(M*.025-1);
	u1[0] = postsamples1(M*.975-1);

	acc(0) = acc(0)/M;
	result = Rcpp::List::create(Rcpp::Named("l0") = l0,Rcpp::Named("u0") = u0,Rcpp::Named("l1") = l1,Rcpp::Named("u1") = u1, Rcpp::Named("acceptance_rate") = acc, Rcpp::Named("samples0") = postsamples0, Rcpp::Named("samples1") = postsamples1, Rcpp::Named("logpost") = logpost);

	return result;
}


inline std::vector<double> GibbsMCMCpsmooth(RVector<double> nn, RMatrix<double> data1, RMatrix<double> data2, RMatrix<double> thetaboot,
	RVector<double> bootmean0, RVector<double> bootmean1, RMatrix<double> databoot1, RMatrix<double> databoot2, RMatrix<double> priordata1, RMatrix<double> priordata2, 
	RVector<double> priorweight, RVector<double> scheduleLen, RMatrix<double> propSched, RVector<double> ddelta, RVector<double> alpha, RVector<double> M_samp, 
	RVector<double> w, std::size_t i) {
   	
	double pi = 3.14159265358979323846;
	int M = int(M_samp[0]);
	int n1 = int(nn[0]);
	int n2 = int(nn[1]);
	int n1p = int(nn[2]);
	int n2p = int(nn[3]);
	int sN = int(scheduleLen[0]);
	double delta0 = double(ddelta[0]);
	double delta1 = double(ddelta[1]);
	double w0 = double(w[0]);
	double w1 = double(w[1]);
   	NumericVector prop0(1,0.0);
   	NumericVector prop1(1,0.0);
	std::vector<double> cov_ind;
	double z1 = 0.0; double z2 = 0.0;
   	NumericVector theta0old(1,0.0);
	NumericVector theta0new(1,0.0);
	NumericVector theta1old(1,0.0);
	NumericVector theta1new(1,0.0);
	NumericVector theta0temp(1,0.0);
	NumericVector theta1temp(1,0.0);
	NumericVector loglikdiff(1,0.0);
	NumericVector loss1temp(1,0.0);
	NumericVector loss2temp(1,0.0);	
	NumericVector loss1old(1,0.0);
	NumericVector loss1new(1,0.0);
	NumericVector loss2old(1,0.0);
	NumericVector loss2new(1,0.0);
	NumericVector loss1tempp(1,0.0);
	NumericVector loss2tempp(1,0.0);	
	NumericVector loss1oldp(1,0.0);
	NumericVector loss1newp(1,0.0);
	NumericVector loss2oldp(1,0.0);
	NumericVector loss2newp(1,0.0);
	NumericVector r(1,0.0);
	NumericVector uu(1,0.0);
	NumericVector postsamples0(M,0.0);
	NumericVector postsamples1(M,0.0);
	NumericVector l0(1,0.0);
	NumericVector l1(1,0.0);
	NumericVector u0(1,0.0);
	NumericVector u1(1,0.0);
	theta0old(0) = thetaboot(i,0);
	theta1old(0) = thetaboot(i,1);
	NumericVector d1y1(1,0.0);
	NumericVector d1y2(1,0.0);
	NumericVector d2y1(1,0.0);
	NumericVector d2y2(1,0.0);
	NumericVector d1y1p(1,0.0);
	NumericVector d1y2p(1,0.0);
	NumericVector d2y1p(1,0.0);
	NumericVector d2y2p(1,0.0);



	
	for(int j=0; j<sN; j++){
		if(w0<=propSched(j,0)){
			prop0(0) = propSched(j,1);	
		}
		if(w1<=propSched(j,2)){
			prop1(0) = propSched(j,3);	
		}		
	}
	if(w0>1.0){
		prop0(0) = propSched(0,1);	
	}
	if(w1>1.0){
		prop1(0) = propSched(0,3);	
	}
	for(int k=0; k<n1; k++){
		if(databoot1(k,2*i)==-1){
			d1y1(0) = d1y1(0) + 1.0;
		}else {
			d1y2(0) = d1y2(0) + 1.0;	
		}
	}
	for(int k=0; k<n2; k++){
		if(databoot2(k,2*i)==-1){
			d2y1(0) = d2y1(0) + 1.0;
		}else {
			d2y2(0) = d2y2(0) + 1.0;	
		}
	}
	for(int k=0; k<n1p; k++){
		if(priordata1(k,0)==-1){
			d1y1p(0) = d1y1p(0) + 1.0;
		}else {
			d1y2p(0) = d1y2p(0) + 1.0;	
		}
	}
	for(int k=0; k<n2p; k++){
		if(priordata2(k,0)==-1){
			d2y1p(0) = d2y1p(0) + 1.0;
		}else {
			d2y2p(0) = d2y2p(0) + 1.0;	
		}
	}

	for(int k=0; k<n1; k++){
		z1 = databoot1(k,2*i)*(databoot1(k,2*i+1)-theta0old(0));
		loss1temp(0)=0;
		if((z1>-delta0/2) && (z1<=delta0/2)){
			loss1temp(0) = 0.5 + 0.5*sin((pi/delta0)*(z1)+pi);
		}else if(z1<0){
			loss1temp(0) = 1.0;	
		}
		if(databoot1(k,2*i)==-1){
			loss1old(0) = loss1old(0) + loss1temp(0)/d1y1(0);
		}else {
			loss1old(0) = loss1old(0) + loss1temp(0)/d1y2(0);
		}
	}
	for(int k=0; k<n2; k++){		
		z2 = databoot2(k,2*i)*(databoot2(k,2*i+1)-theta1old(0));
		loss2temp(0)=0;
		if((z2>-delta1/2) && (z2<=delta1/2)){
			loss2temp(0) = 0.5 + 0.5*sin((pi/delta1)*(z2)+pi);
		}else if(z2<0){
			loss2temp(0) = 1.0;	
		}
		if(databoot2(k,2*i)==-1){
			loss2old(0) = loss2old(0) + loss2temp(0)/d2y1(0);
		}else {
			loss2old(0) = loss2old(0) + loss2temp(0)/d2y2(0);
		}
	}
	
	for(int k=0; k<n1p; k++){
		z1 = priordata1(k,0)*(priordata1(k,1)-theta0old(0));
		loss1tempp(0)=0;
		if((z1>-delta0/2) && (z1<=delta0/2)){
			loss1tempp(0) = 0.5 + 0.5*sin((pi/delta0)*(z1)+pi);
		}else if(z1<0){
			loss1tempp(0) = 1.0;	
		}
		if(priordata1(k,0)==-1){
			loss1oldp(0) = loss1oldp(0) + loss1tempp(0)/d1y1p(0);
		}else {
			loss1oldp(0) = loss1oldp(0) + loss1tempp(0)/d1y2p(0);
		}
	}
	for(int k=0; k<n2p; k++){		
		z2 = priordata2(k,0)*(priordata2(k,1)-theta1old(0));
		loss2tempp(0)=0;
		if((z2>-delta1/2) && (z2<=delta1/2)){
			loss2tempp(0) = 0.5 + 0.5*sin((pi/delta1)*(z2)+pi);
		}else if(z2<0){
			loss2tempp(0) = 1.0;	
		}
		if(priordata2(k,0)==-1){
			loss2oldp(0) = loss2oldp(0) + loss2tempp(0)/d2y1p(0);
		}else {
			loss2oldp(0) = loss2oldp(0) + loss2tempp(0)/d2y2p(0);
		}
	}

	for(int j=0; j<M; j++) {
		theta0temp(0) = R::rnorm(theta0old(0), prop0(0));
		theta1temp(0) = R::rnorm(theta1old(0), prop1(0));
		theta0new(0) = min(theta0temp(0), theta1temp(0));
		theta1new(0) = max(theta0temp(0), theta1temp(0));
		
		loss1new(0) = 0.0;
		loss2new(0) = 0.0;
		loss1newp(0) = 0.0;
		loss2newp(0) = 0.0;
		loglikdiff(0) = 0.0;
		for(int k=0; k<n1; k++){
			z1 = databoot1(k,2*i)*(databoot1(k,2*i+1)-theta0new(0));
			loss1temp(0)=0;
			if((z1>-delta0/2) && (z1<=delta0/2)){
				loss1temp(0) = 0.5 + 0.5*sin((pi/delta0)*(z1)+pi);
			}else if(z1<0){
				loss1temp(0) = 1.0;	
			}
			if(databoot1(k,2*i)==-1){
				loss1new(0) = loss1new(0) + loss1temp(0)/d1y1(0);
			}else {
				loss1new(0) = loss1new(0) + loss1temp(0)/d1y2(0);
			}
		}
		for(int k=0; k<n2; k++){		
			z2 = databoot2(k,2*i)*(databoot2(k,2*i+1)-theta1new(0));
			loss2temp(0)=0;
			if((z2>-delta1/2) && (z2<=delta1/2)){
				loss2temp(0) = 0.5 + 0.5*sin((pi/delta1)*(z2)+pi);
			}else if(z2<0){
				loss2temp(0) = 1.0;	
			}
			if(databoot2(k,2*i)==-1){
				loss2new(0) = loss2new(0) + loss2temp(0)/d2y1(0);
			}else {
				loss2new(0) = loss2new(0) + loss2temp(0)/d2y2(0);
			}
		}
		
		for(int k=0; k<n1p; k++){
			z1 = priordata1(k,0)*(priordata1(k,1)-theta0new(0));
			loss1tempp(0)=0;
			if((z1>-delta0/2) && (z1<=delta0/2)){
				loss1tempp(0) = 0.5 + 0.5*sin((pi/delta0)*(z1)+pi);
			}else if(z1<0){
				loss1tempp(0) = 1.0;	
			}
			if(priordata1(k,0)==-1){
				loss1newp(0) = loss1newp(0) + loss1tempp(0)/d1y1p(0);
			}else {
				loss1newp(0) = loss1newp(0) + loss1tempp(0)/d1y2p(0);
			}
		}
		for(int k=0; k<n2p; k++){		
			z2 = priordata2(k,0)*(priordata2(k,1)-theta1new(0));
			loss2tempp(0)=0;
			if((z2>-delta1/2) && (z2<=delta1/2)){
				loss2tempp(0) = 0.5 + 0.5*sin((pi/delta1)*(z2)+pi);
			}else if(z2<0){
				loss2tempp(0) = 1.0;	
			}
			if(priordata2(k,0)==-1){
				loss2newp(0) = loss2newp(0) + loss2tempp(0)/d2y1p(0);
			}else {
				loss2newp(0) = loss2newp(0) + loss2tempp(0)/d2y2p(0);
			}
		}

		loglikdiff(0) = -n2*w1*(loss2new(0)-loss2old(0))-n1*w0*(loss1new(0)-loss1old(0))-n2p*priorweight[1]*(loss2newp(0)-loss2oldp(0))-n1p*priorweight[0]*(loss1newp(0)-loss1oldp(0));
		loglikdiff(0) = fmin(std::exp(loglikdiff(0)), 1.0);
		uu[0] = R::runif(0.0,1.0);
		if(uu(0) <= loglikdiff(0)) {
			postsamples0(j) = theta0new(0);
			postsamples1(j) = theta1new(0);
			theta0old(0) = theta0new(0);
			theta1old(0) = theta1new(0);
			loss2old(0) = loss2new(0);
			loss1old(0) = loss1new(0);	
		}
		else {
			postsamples0(j) = theta0old(0);
			postsamples1(j) = theta1old(0);
		}
	}	
/*
	double templogpost;
	double tempYI;
	double temppost0;
	double temppost1;
	bool swapped;	
        for (int j = 0; j < M-1; j++){ 
		 swapped = false; 
       		 for (int k = 0; k < M-j-1; k++){  
           		 if (logpost(k) > logpost(k+1)){ 
				 swapped = true;
				 templogpost = logpost(k);
				 logpost(k) = logpost(k+1);
				 logpost(k+1) = templogpost;
				 tempYI = YI(k);
				 YI(k) = YI(k+1);
				 YI(k+1) = tempYI;
				 temppost0 = postsamples0(k);
				 postsamples0(k) = postsamples0(k+1);
				 postsamples0(k+1) = temppost0;				 
				 temppost1 = postsamples1(k);
				 postsamples1(k) = postsamples1(k+1);
				 postsamples1(k+1) = temppost1;
			 }
		 }
		 if (swapped == false) 
			break;
        }
			
	l0[0] = 100000;
	u0[0] = -100000;
	l1[0] = 100000;
	u1[0] = -100000;
	YIl[0] = 100000;
	YIu[0] = -100000;
	for(int i = (0.05*M-1); i<M; i++){
		if(l0[0]>postsamples0(i)){
			l0[0]=postsamples0(i);
		}
		if(u0[0]<postsamples0(i)){
			u0[0]=postsamples0(i);
		}
		if(l1[0]>postsamples1(i)){
			l1[0]=postsamples1(i);
		}
		if(u1[0]<postsamples1(i)){
			u1[0]=postsamples1(i);
		}
		if(YIl[0]>YI(i)){
			YIl[0]=YI(i);
		}
		if(YIu[0]<YI(i)){
			YIu[0]=YI(i);
		}
	}*/
	
	std::sort(postsamples0.begin(), postsamples0.end());
	std::sort(postsamples1.begin(), postsamples1.end());
	l0[0] = postsamples0(M*.025-1);
	u0[0] = postsamples0(M*.975-1);
	l1[0] = postsamples1(M*.025-1);
	u1[0] = postsamples1(M*.975-1);
	if ( (l0[0] < bootmean0[0]) && (u0[0] > bootmean0[0]) ){
		cov_ind.push_back(1.0);
	} else {cov_ind.push_back(0.0);}
	if ( (l1[0] < bootmean1[0]) && (u1[0] > bootmean1[0]) ){
		cov_ind.push_back(1.0);
	} else {cov_ind.push_back(0.0);}
	
	return cov_ind;	
}


// Proposal Schedule as input
// [[Rcpp::export]]
Rcpp::List GibbsMCMCp2smooth(NumericVector nn, NumericMatrix data1, NumericMatrix data2, NumericVector bootmean0, NumericVector bootmean1, NumericMatrix priordata1, NumericMatrix priordata2, NumericVector priorweight, NumericVector scheduleLen, NumericMatrix propSched, NumericVector alpha, NumericVector ddelta, NumericVector M_samp, NumericVector w) {
	
	List result;
	double pi = 3.14159265358979323846;
	int M = int(M_samp[0]);
	int n1 = int(nn[0]);
	int n2 = int(nn[1]);
	int n1p = int(nn[2]);
	int n2p = int(nn[3]);
	int sN = int(scheduleLen[0]);
	double delta0 = double(ddelta[0]);
	double delta1 = double(ddelta[1]);
	double w0 = double(w[0]);
	double w1 = double(w[1]);
	double pw1 = double(priorweight[0]);
	double pw2 = double(priorweight[1]);
   	NumericVector prop0(1,0.0);
   	NumericVector prop1(1,0.0);
	double z1 = 0.0; double z2 = 0.0;
   	NumericVector theta0old(1,0.0);
	NumericVector theta0new(1,0.0);
	NumericVector theta1old(1,0.0);
	NumericVector theta1new(1,0.0);
	NumericVector theta0temp(1,0.0);
	NumericVector theta1temp(1,0.0);
	NumericVector loglikdiff(1,0.0);
	NumericVector loglikdiff1(1,0.0);
	NumericVector loss1temp(1,0.0);
	NumericVector loss2temp(1,0.0);	
	NumericVector loss1old(1,0.0);
	NumericVector loss1new(1,0.0);
	NumericVector loss2old(1,0.0);
	NumericVector loss2new(1,0.0);
	NumericVector loss1tempp(1,0.0);
	NumericVector loss2tempp(1,0.0);	
	NumericVector loss1oldp(1,0.0);
	NumericVector loss1newp(1,0.0);
	NumericVector loss2oldp(1,0.0);
	NumericVector loss2newp(1,0.0);
	NumericVector r(1,0.0);
	NumericVector uu(1,0.0);
	NumericVector postsamples0(M,0.0);
	NumericVector postsamples1(M,0.0);
	NumericVector logpost(M,0.0);
	NumericVector l0(1,0.0);
	NumericVector l1(1,0.0);
	NumericVector u0(1,0.0);
	NumericVector u1(1,0.0);
	theta0old(0) = bootmean0(0);
	theta1old(0) = bootmean1(0);
	NumericVector acc(1, 0.0);
	NumericVector d1y1(1,0.0);
	NumericVector d1y2(1,0.0);
	NumericVector d2y1(1,0.0);
	NumericVector d2y2(1,0.0);
	NumericVector d1y1p(1,0.0);
	NumericVector d1y2p(1,0.0);
	NumericVector d2y1p(1,0.0);
	NumericVector d2y2p(1,0.0);


	
	for(int j=0; j<sN; j++){
		if(w0<=propSched(j,0)){
			prop0(0) = propSched(j,1);	
		}
		if(w1<=propSched(j,2)){
			prop1(0) = propSched(j,3);	
		}		
	}
	if(w0>1.0){
		prop0(0) = propSched(0,1);	
	}
	if(w1>1.0){
		prop1(0) = propSched(0,3);	
	}
	
	for(int k=0; k<n1; k++){
		if(data1(k,0)==-1){
			d1y1(0) = d1y1(0) + 1.0;
		}else {
			d1y2(0) = d1y2(0) + 1.0;	
		}
	}
	for(int k=0; k<n2; k++){
		if(data2(k,0)==-1){
			d2y1(0) = d2y1(0) + 1.0;
		}else {
			d2y2(0) = d2y2(0) + 1.0;	
		}
	}
	for(int k=0; k<n1p; k++){
		if(priordata1(k,0)==-1){
			d1y1p(0) = d1y1p(0) + 1.0;
		}else {
			d1y2p(0) = d1y2p(0) + 1.0;	
		}
	}
	for(int k=0; k<n2p; k++){
		if(priordata2(k,0)==-1){
			d2y1p(0) = d2y1p(0) + 1.0;
		}else {
			d2y2p(0) = d2y2p(0) + 1.0;	
		}
	}
	
	for(int k=0; k<n1; k++){
		z1 = data1(k,0)*(data1(k,1)-theta0old(0));
		loss1temp(0)=0;
		if((z1>-delta0/2) && (z1<=delta0/2)){
			loss1temp(0) = 0.5 + 0.5*sin((pi/delta0)*(z1)+pi);
		}else if(z1<0){
			loss1temp(0) = 1.0;	
		}
		if(data1(k,0)==-1){
			loss1old(0) = loss1old(0) + loss1temp(0)/d1y1(0);
		}else {
			loss1old(0) = loss1old(0) + loss1temp(0)/d1y2(0);
		}
	}
	for(int k=0; k<n2; k++){
		z2 = data2(k,0)*(data2(k,1)-theta1old(0));
		loss2temp(0)=0;
		if((z2>-delta1/2) && (z2<=delta1/2)){
			loss2temp(0) = 0.5 + 0.5*sin((pi/delta1)*(z2)+pi);
		}else if(z2<0){
			loss2temp(0) = 1.0;	
		}
		if(data2(k,0)==-1){
			loss2old(0) = loss2old(0) + loss2temp(0)/d2y1(0);
		}else {
			loss2old(0) = loss2old(0) + loss2temp(0)/d2y2(0);
		}
	}


	for(int k=0; k<n1p; k++){
		z1 = priordata1(k,0)*(priordata1(k,1)-theta0old(0));
		loss1tempp(0)=0;
		if((z1>-delta0/2) && (z1<=delta0/2)){
			loss1tempp(0) = 0.5 + 0.5*sin((pi/delta0)*(z1)+pi);
		}else if(z1<0){
			loss1tempp(0) = 1.0;	
		}
		if(priordata1(k,0)==-1){
			loss1oldp(0) = loss1oldp(0) + loss1tempp(0)/d1y1p(0);
		}else {
			loss1oldp(0) = loss1oldp(0) + loss1tempp(0)/d1y2p(0);
		}
	}
	for(int k=0; k<n2p; k++){
		z2 = priordata2(k,0)*(priordata2(k,1)-theta1old(0));
		loss2tempp(0)=0;
		if((z2>-delta1/2) && (z2<=delta1/2)){
			loss2tempp(0) = 0.5 + 0.5*sin((pi/delta1)*(z2)+pi);
		}else if(z2<0){
			loss2tempp(0) = 1.0;	
		}
		if(priordata2(k,0)==-1){
			loss2oldp(0) = loss2oldp(0) + loss2tempp(0)/d2y1p(0);
		}else {
			loss2oldp(0) = loss2oldp(0) + loss2tempp(0)/d2y2p(0);
		}
	}

	for(int j=0; j<M; j++) {
		theta0temp(0) = R::rnorm(theta0old(0), prop0(0));
		theta1temp(0) = R::rnorm(theta1old(0), prop1(0));
		theta0new(0) = min(theta0temp(0), theta1temp(0));
		theta1new(0) = max(theta0temp(0), theta1temp(0));

		
		loss1new(0) = 0.0;
		loss2new(0) = 0.0;
		loss1newp(0) = 0.0;
		loss2newp(0) = 0.0;
		for(int k=0; k<n1; k++){
			z1 = data1(k,0)*(data1(k,1)-theta0new(0));
			loss1temp(0)=0;
			if((z1>-delta0/2) && (z1<=delta0/2)){
				loss1temp(0) = 0.5 + 0.5*sin((pi/delta0)*(z1)+pi);
			}else if(z1<0){
				loss1temp(0) = 1.0;	
			}
			if(data1(k,0)==-1){
				loss1new(0) = loss1new(0) + loss1temp(0)/d1y1(0);
			}else {
				loss1new(0) = loss1new(0) + loss1temp(0)/d1y2(0);
			}
		}
		for(int k=0; k<n2; k++){
			z2 = data2(k,0)*(data2(k,1)-theta1new(0));
			loss2temp(0)=0;
			if((z2>-delta1/2) && (z2<=delta1/2)){
				loss2temp(0) = 0.5 + 0.5*sin((pi/delta1)*(z2)+pi);
			}else if(z2<0){
				loss2temp(0) = 1.0;	
			}
			if(data2(k,0)==-1){
				loss2new(0) = loss2new(0) + loss2temp(0)/d2y1(0);
			}else {
				loss2new(0) = loss2new(0) + loss2temp(0)/d2y2(0);
			}
		}
		for(int k=0; k<n1p; k++){
			z1 = priordata1(k,0)*(priordata1(k,1)-theta0new(0));
			loss1tempp(0)=0;
			if((z1>-delta0/2) && (z1<=delta0/2)){
				loss1tempp(0) = 0.5 + 0.5*sin((pi/delta0)*(z1)+pi);
			}else if(z1<0){
				loss1tempp(0) = 1.0;	
			}
			if(priordata1(k,0)==-1){
				loss1newp(0) = loss1newp(0) + loss1tempp(0)/d1y1p(0);
			}else {
				loss1newp(0) = loss1newp(0) + loss1tempp(0)/d1y2p(0);
			}
		}
		for(int k=0; k<n2p; k++){
			z2 = priordata2(k,0)*(priordata2(k,1)-theta1new(0));
			loss2tempp(0)=0;
			if((z2>-delta1/2) && (z2<=delta1/2)){
				loss2tempp(0) = 0.5 + 0.5*sin((pi/delta1)*(z2)+pi);
			}else if(z2<0){
				loss2tempp(0) = 1.0;	
			}
			if(priordata2(k,0)==-1){
				loss2newp(0) = loss2newp(0) + loss2tempp(0)/d2y1p(0);
			}else {
				loss2newp(0) = loss2newp(0) + loss2tempp(0)/d2y2p(0);
			}
		}
		loglikdiff(0) = 0.0;
		loglikdiff1(0) = 0.0;
		loglikdiff(0) = -n2*w1*(loss2new(0)-loss2old(0))-n1*w0*(loss1new(0)-loss1old(0)) -n2p*pw2*(loss2newp(0)-loss2oldp(0))-n1p*pw1*(loss1newp(0)-loss1oldp(0));
		loglikdiff1(0) = fmin(std::exp(loglikdiff(0)), 1.0);
		uu[0] = R::runif(0.0,1.0);
		if(uu(0) <= loglikdiff1(0)) {
			postsamples0(j) = theta0new(0);
			postsamples1(j) = theta1new(0);
			logpost(j) = log(loglikdiff(0));
			theta0old(0) = theta0new(0);
			theta1old(0) = theta1new(0);
			loss2old(0) = loss2new(0);
			loss1old(0) = loss1new(0);
			loss2oldp(0) = loss2newp(0);
			loss1oldp(0) = loss1newp(0);
			acc(0) = acc(0)+1;
		}
		else {
			postsamples0(j) = theta0old(0);
			postsamples1(j) = theta1old(0);
			if(j>0){
				logpost(j) = logpost(j-1);
			}else {
				logpost(j) = logpost(0);	
			}
		}
	}	

	std::sort(postsamples0.begin(), postsamples0.end());
	std::sort(postsamples1.begin(), postsamples1.end());
	l0[0] = postsamples0(M*.025-1);
	u0[0] = postsamples0(M*.975-1);
	l1[0] = postsamples1(M*.025-1);
	u1[0] = postsamples1(M*.975-1);

	acc(0) = acc(0)/M;
	result = Rcpp::List::create(Rcpp::Named("l0") = l0,Rcpp::Named("u0") = u0,Rcpp::Named("l1") = l1,Rcpp::Named("u1") = u1, Rcpp::Named("acceptance_rate") = acc, Rcpp::Named("samples0") = postsamples0, Rcpp::Named("samples1") = postsamples1, Rcpp::Named("logpost") = logpost);

	
	return result;
}


inline std::vector<double> GibbsMCMC3loopsmooth(RVector<double> nn, RMatrix<double> data1, RMatrix<double> data2, RMatrix<double> thetaboot,
	RVector<double> bootmean0, RVector<double> bootmean1, RMatrix<double> databoot1, RMatrix<double> databoot2, 
	RVector<double> priorsched, RVector<double> scheduleLen, RMatrix<double> propSched, RVector<double> ddelta, RVector<double> alpha, RVector<double> M_samp, 
	RVector<double> w, std::size_t i) {
   	
	double pi = 3.14159265358979323846;
	int M = int(M_samp[0]);
	int n1 = int(nn[0]);
	int n2 = int(nn[1]);
	int sN = int(scheduleLen[0]);
	double delta0 = double(ddelta[0]);
	double delta1 = double(ddelta[1]);
	double w0 = double(w[0]);
	double w1 = double(w[1]);
   	NumericVector prop0(1,0.0);
   	NumericVector prop1(1,0.0);
	std::vector<double> cov_ind;
	double z1 = 0.0; double z2 = 0.0;
   	NumericVector theta0old(1,0.0);
	NumericVector theta0new(1,0.0);
	NumericVector theta1old(1,0.0);
	NumericVector theta1new(1,0.0);
	NumericVector theta0temp(1,0.0);
	NumericVector theta1temp(1,0.0);
	NumericVector loglikdiff(1,0.0);
	NumericVector loss1temp(1,0.0);
	NumericVector loss2temp(1,0.0);	
	NumericVector loss1old(1,0.0);
	NumericVector loss1new(1,0.0);
	NumericVector loss2old(1,0.0);
	NumericVector loss2new(1,0.0);
	NumericVector r(1,0.0);
	NumericVector uu(1,0.0);
	NumericVector postsamples0(M,0.0);
	NumericVector postsamples1(M,0.0);
	NumericVector l0(1,0.0);
	NumericVector l1(1,0.0);
	NumericVector u0(1,0.0);
	NumericVector u1(1,0.0);
	theta0old(0) = thetaboot(i,0);
	theta1old(0) = thetaboot(i,1);
	NumericVector d1y1(1,0.0);
	NumericVector d1y2(1,0.0);
	NumericVector d2y1(1,0.0);
	NumericVector d2y2(1,0.0);




	
	for(int j=0; j<sN; j++){
		if(w0<=propSched(j,0)){
			prop0(0) = propSched(j,1);	
		}
		if(w1<=propSched(j,2)){
			prop1(0) = propSched(j,3);	
		}		
	}
	if(w0>1.0){
		prop0(0) = propSched(0,1);	
	}
	if(w1>1.0){
		prop1(0) = propSched(0,3);	
	}
	for(int k=0; k<n1; k++){
		if(databoot1(k,2*i)==-1){
			d1y1(0) = d1y1(0) + 1.0;
		}else {
			d1y2(0) = d1y2(0) + 1.0;	
		}
	}
	for(int k=0; k<n2; k++){
		if(databoot2(k,2*i)==-1){
			d2y1(0) = d2y1(0) + 1.0;
		}else {
			d2y2(0) = d2y2(0) + 1.0;	
		}
	}
	

	for(int k=0; k<n1; k++){
		z1 = databoot1(k,2*i)*(databoot1(k,2*i+1)-theta0old(0));
		loss1temp(0)=0;
		if((z1>-delta0/2) && (z1<=delta0/2)){
			loss1temp(0) = 0.5 + 0.5*sin((pi/delta0)*(z1)+pi);
		}else if(z1<0){
			loss1temp(0) = 1.0;	
		}
		if(databoot1(k,2*i)==-1){
			loss1old(0) = loss1old(0) + loss1temp(0)/d1y1(0);
		}else {
			loss1old(0) = loss1old(0) + loss1temp(0)/d1y2(0);
		}
	}
	for(int k=0; k<n2; k++){		
		z2 = databoot2(k,2*i)*(databoot2(k,2*i+1)-theta1old(0));
		loss2temp(0)=0;
		if((z2>-delta1/2) && (z2<=delta1/2)){
			loss2temp(0) = 0.5 + 0.5*sin((pi/delta1)*(z2)+pi);
		}else if(z2<0){
			loss2temp(0) = 1.0;	
		}
		if(databoot2(k,2*i)==-1){
			loss2old(0) = loss2old(0) + loss2temp(0)/d2y1(0);
		}else {
			loss2old(0) = loss2old(0) + loss2temp(0)/d2y2(0);
		}
	}
	
	

	for(int j=0; j<M; j++) {
		theta0temp(0) = R::rnorm(theta0old(0), prop0(0));
		theta1temp(0) = R::rnorm(theta1old(0), prop1(0));
		theta0new(0) = min(theta0temp(0), theta1temp(0));
		theta1new(0) = max(theta0temp(0), theta1temp(0));
		
		loss1new(0) = 0.0;
		loss2new(0) = 0.0;
		
		loglikdiff(0) = 0.0;
		for(int k=0; k<n1; k++){
			z1 = databoot1(k,2*i)*(databoot1(k,2*i+1)-theta0new(0));
			loss1temp(0)=0;
			if((z1>-delta0/2) && (z1<=delta0/2)){
				loss1temp(0) = 0.5 + 0.5*sin((pi/delta0)*(z1)+pi);
			}else if(z1<0){
				loss1temp(0) = 1.0;	
			}
			if(databoot1(k,2*i)==-1){
				loss1new(0) = loss1new(0) + loss1temp(0)/d1y1(0);
			}else {
				loss1new(0) = loss1new(0) + loss1temp(0)/d1y2(0);
			}
		}
		for(int k=0; k<n2; k++){		
			z2 = databoot2(k,2*i)*(databoot2(k,2*i+1)-theta1new(0));
			loss2temp(0)=0;
			if((z2>-delta1/2) && (z2<=delta1/2)){
				loss2temp(0) = 0.5 + 0.5*sin((pi/delta1)*(z2)+pi);
			}else if(z2<0){
				loss2temp(0) = 1.0;	
			}
			if(databoot2(k,2*i)==-1){
				loss2new(0) = loss2new(0) + loss2temp(0)/d2y1(0);
			}else {
				loss2new(0) = loss2new(0) + loss2temp(0)/d2y2(0);
			}
		}
		
		

		loglikdiff(0) = -n2*w1*(loss2new(0)-loss2old(0))-n1*w0*(loss1new(0)-loss1old(0));
		loglikdiff(0) = fmin(std::exp(loglikdiff(0))*((R::dnorm(theta0new(0),priorsched[0],priorsched[1],0)*R::dnorm(theta1new(0),priorsched[2],priorsched[3],0)+R::dnorm(theta1new(0),priorsched[0],priorsched[1],0)*R::dnorm(theta0new(0),priorsched[2],priorsched[3],0))/(R::dnorm(theta0old(0),priorsched[0],priorsched[1],0)*R::dnorm(theta1old(0),priorsched[2],priorsched[3],0)+R::dnorm(theta1old(0),priorsched[0],priorsched[1],0)*R::dnorm(theta0old(0),priorsched[2],priorsched[3],0))), 1.0);
		uu[0] = R::runif(0.0,1.0);
		if(uu(0) <= loglikdiff(0)) {
			postsamples0(j) = theta0new(0);
			postsamples1(j) = theta1new(0);
			theta0old(0) = theta0new(0);
			theta1old(0) = theta1new(0);
			loss2old(0) = loss2new(0);
			loss1old(0) = loss1new(0);	
		}
		else {
			postsamples0(j) = theta0old(0);
			postsamples1(j) = theta1old(0);
		}
	}	
/*
	double templogpost;
	double tempYI;
	double temppost0;
	double temppost1;
	bool swapped;	
        for (int j = 0; j < M-1; j++){ 
		 swapped = false; 
       		 for (int k = 0; k < M-j-1; k++){  
           		 if (logpost(k) > logpost(k+1)){ 
				 swapped = true;
				 templogpost = logpost(k);
				 logpost(k) = logpost(k+1);
				 logpost(k+1) = templogpost;
				 tempYI = YI(k);
				 YI(k) = YI(k+1);
				 YI(k+1) = tempYI;
				 temppost0 = postsamples0(k);
				 postsamples0(k) = postsamples0(k+1);
				 postsamples0(k+1) = temppost0;				 
				 temppost1 = postsamples1(k);
				 postsamples1(k) = postsamples1(k+1);
				 postsamples1(k+1) = temppost1;
			 }
		 }
		 if (swapped == false) 
			break;
        }
			
	l0[0] = 100000;
	u0[0] = -100000;
	l1[0] = 100000;
	u1[0] = -100000;
	YIl[0] = 100000;
	YIu[0] = -100000;
	for(int i = (0.05*M-1); i<M; i++){
		if(l0[0]>postsamples0(i)){
			l0[0]=postsamples0(i);
		}
		if(u0[0]<postsamples0(i)){
			u0[0]=postsamples0(i);
		}
		if(l1[0]>postsamples1(i)){
			l1[0]=postsamples1(i);
		}
		if(u1[0]<postsamples1(i)){
			u1[0]=postsamples1(i);
		}
		if(YIl[0]>YI(i)){
			YIl[0]=YI(i);
		}
		if(YIu[0]<YI(i)){
			YIu[0]=YI(i);
		}
	}*/
	
	std::sort(postsamples0.begin(), postsamples0.end());
	std::sort(postsamples1.begin(), postsamples1.end());
	l0[0] = postsamples0(M*.025-1);
	u0[0] = postsamples0(M*.975-1);
	l1[0] = postsamples1(M*.025-1);
	u1[0] = postsamples1(M*.975-1);
	if ( (l0[0] < bootmean0[0]) && (u0[0] > bootmean0[0]) ){
		cov_ind.push_back(1.0);
	} else {cov_ind.push_back(0.0);}
	if ( (l1[0] < bootmean1[0]) && (u1[0] > bootmean1[0]) ){
		cov_ind.push_back(1.0);
	} else {cov_ind.push_back(0.0);}
	
	return cov_ind;	
}


// Proposal Schedule as input
// [[Rcpp::export]]
Rcpp::List GibbsMCMC3smooth(NumericVector nn, NumericMatrix data1, NumericMatrix data2, NumericVector bootmean0, NumericVector bootmean1, NumericVector priorsched, NumericVector scheduleLen, NumericMatrix propSched, NumericVector alpha, NumericVector ddelta, NumericVector M_samp, NumericVector w) {
	
	List result;
	double pi = 3.14159265358979323846;
	int M = int(M_samp[0]);
	int n1 = int(nn[0]);
	int n2 = int(nn[1]);
	int sN = int(scheduleLen[0]);
	double delta0 = double(ddelta[0]);
	double delta1 = double(ddelta[1]);
	double w0 = double(w[0]);
	double w1 = double(w[1]);
   	NumericVector prop0(1,0.0);
   	NumericVector prop1(1,0.0);
	double z1 = 0.0; double z2 = 0.0;
   	NumericVector theta0old(1,0.0);
	NumericVector theta0new(1,0.0);
	NumericVector theta1old(1,0.0);
	NumericVector theta1new(1,0.0);
	NumericVector theta0temp(1,0.0);
	NumericVector theta1temp(1,0.0);
	NumericVector loglikdiff(1,0.0);
	NumericVector loglikdiff1(1,0.0);
	NumericVector loss1temp(1,0.0);
	NumericVector loss2temp(1,0.0);	
	NumericVector loss1old(1,0.0);
	NumericVector loss1new(1,0.0);
	NumericVector loss2old(1,0.0);
	NumericVector loss2new(1,0.0);
	NumericVector r(1,0.0);
	NumericVector uu(1,0.0);
	NumericVector postsamples0(M,0.0);
	NumericVector postsamples1(M,0.0);
	NumericVector logpost(M,0.0);
	NumericVector l0(1,0.0);
	NumericVector l1(1,0.0);
	NumericVector u0(1,0.0);
	NumericVector u1(1,0.0);
	theta0old(0) = bootmean0(0);
	theta1old(0) = bootmean1(0);
	NumericVector acc(1, 0.0);
	NumericVector d1y1(1,0.0);
	NumericVector d1y2(1,0.0);
	NumericVector d2y1(1,0.0);
	NumericVector d2y2(1,0.0);



	
	for(int j=0; j<sN; j++){
		if(w0<=propSched(j,0)){
			prop0(0) = propSched(j,1);	
		}
		if(w1<=propSched(j,2)){
			prop1(0) = propSched(j,3);	
		}		
	}
	if(w0>1.0){
		prop0(0) = propSched(0,1);	
	}
	if(w1>1.0){
		prop1(0) = propSched(0,3);	
	}
	
	for(int k=0; k<n1; k++){
		if(data1(k,0)==-1){
			d1y1(0) = d1y1(0) + 1.0;
		}else {
			d1y2(0) = d1y2(0) + 1.0;	
		}
	}
	for(int k=0; k<n2; k++){
		if(data2(k,0)==-1){
			d2y1(0) = d2y1(0) + 1.0;
		}else {
			d2y2(0) = d2y2(0) + 1.0;	
		}
	}

	
	for(int k=0; k<n1; k++){
		z1 = data1(k,0)*(data1(k,1)-theta0old(0));
		loss1temp(0)=0;
		if((z1>-delta0/2) && (z1<=delta0/2)){
			loss1temp(0) = 0.5 + 0.5*sin((pi/delta0)*(z1)+pi);
		}else if(z1<0){
			loss1temp(0) = 1.0;	
		}
		if(data1(k,0)==-1){
			loss1old(0) = loss1old(0) + loss1temp(0)/d1y1(0);
		}else {
			loss1old(0) = loss1old(0) + loss1temp(0)/d1y2(0);
		}
	}
	for(int k=0; k<n2; k++){
		z2 = data2(k,0)*(data2(k,1)-theta1old(0));
		loss2temp(0)=0;
		if((z2>-delta1/2) && (z2<=delta1/2)){
			loss2temp(0) = 0.5 + 0.5*sin((pi/delta1)*(z2)+pi);
		}else if(z2<0){
			loss2temp(0) = 1.0;	
		}
		if(data2(k,0)==-1){
			loss2old(0) = loss2old(0) + loss2temp(0)/d2y1(0);
		}else {
			loss2old(0) = loss2old(0) + loss2temp(0)/d2y2(0);
		}
	}



	for(int j=0; j<M; j++) {
		theta0temp(0) = R::rnorm(theta0old(0), prop0(0));
		theta1temp(0) = R::rnorm(theta1old(0), prop1(0));
		theta0new(0) = min(theta0temp(0), theta1temp(0));
		theta1new(0) = max(theta0temp(0), theta1temp(0));

		
		loss1new(0) = 0.0;
		loss2new(0) = 0.0;

		for(int k=0; k<n1; k++){
			z1 = data1(k,0)*(data1(k,1)-theta0new(0));
			loss1temp(0)=0;
			if((z1>-delta0/2) && (z1<=delta0/2)){
				loss1temp(0) = 0.5 + 0.5*sin((pi/delta0)*(z1)+pi);
			}else if(z1<0){
				loss1temp(0) = 1.0;	
			}
			if(data1(k,0)==-1){
				loss1new(0) = loss1new(0) + loss1temp(0)/d1y1(0);
			}else {
				loss1new(0) = loss1new(0) + loss1temp(0)/d1y2(0);
			}
		}
		for(int k=0; k<n2; k++){
			z2 = data2(k,0)*(data2(k,1)-theta1new(0));
			loss2temp(0)=0;
			if((z2>-delta1/2) && (z2<=delta1/2)){
				loss2temp(0) = 0.5 + 0.5*sin((pi/delta1)*(z2)+pi);
			}else if(z2<0){
				loss2temp(0) = 1.0;	
			}
			if(data2(k,0)==-1){
				loss2new(0) = loss2new(0) + loss2temp(0)/d2y1(0);
			}else {
				loss2new(0) = loss2new(0) + loss2temp(0)/d2y2(0);
			}
		}

		loglikdiff(0) = 0.0;
		loglikdiff1(0) = 0.0;
		loglikdiff(0) = -n2*w1*(loss2new(0)-loss2old(0))-n1*w0*(loss1new(0)-loss1old(0));
		loglikdiff1(0) = fmin(std::exp(loglikdiff(0))*((R::dnorm(theta0new(0),priorsched[0],priorsched[1],0)*R::dnorm(theta1new(0),priorsched[2],priorsched[3],0)+R::dnorm(theta1new(0),priorsched[0],priorsched[1],0)*R::dnorm(theta0new(0),priorsched[2],priorsched[3],0))/(R::dnorm(theta0old(0),priorsched[0],priorsched[1],0)*R::dnorm(theta1old(0),priorsched[2],priorsched[3],0)+R::dnorm(theta1old(0),priorsched[0],priorsched[1],0)*R::dnorm(theta0old(0),priorsched[2],priorsched[3],0))), 1.0);
		uu[0] = R::runif(0.0,1.0);
		if(uu(0) <= loglikdiff1(0)) {
			postsamples0(j) = theta0new(0);
			postsamples1(j) = theta1new(0);
			logpost(j) = log(loglikdiff(0));
			theta0old(0) = theta0new(0);
			theta1old(0) = theta1new(0);
			loss2old(0) = loss2new(0);
			loss1old(0) = loss1new(0);

			acc(0) = acc(0)+1;
		}
		else {
			postsamples0(j) = theta0old(0);
			postsamples1(j) = theta1old(0);
			if(j>0){
				logpost(j) = logpost(j-1);
			}else {
				logpost(j) = logpost(0);	
			}
		}
	}	

	std::sort(postsamples0.begin(), postsamples0.end());
	std::sort(postsamples1.begin(), postsamples1.end());
	l0[0] = postsamples0(M*.025-1);
	u0[0] = postsamples0(M*.975-1);
	l1[0] = postsamples1(M*.025-1);
	u1[0] = postsamples1(M*.975-1);

	acc(0) = acc(0)/M;
	result = Rcpp::List::create(Rcpp::Named("l0") = l0,Rcpp::Named("u0") = u0,Rcpp::Named("l1") = l1,Rcpp::Named("u1") = u1, Rcpp::Named("acceptance_rate") = acc, Rcpp::Named("samples0") = postsamples0, Rcpp::Named("samples1") = postsamples1, Rcpp::Named("logpost") = logpost);

	
	return result;
}



// Proposal Schedule as input
// [[Rcpp::export]]
Rcpp::List GibbsMCMC2(NumericVector nn, NumericMatrix data, NumericVector nnp, NumericMatrix priordata, NumericVector priorweight, NumericMatrix thetaboot,
	NumericVector bootmean0, NumericVector bootmean1, NumericVector scheduleLen, NumericMatrix priorSched, NumericVector alpha, NumericVector M_samp, NumericVector w) {
   	
	List result;
	int M = int(M_samp[0]);
	int n = int(nn[0]);
	int sN = int(scheduleLen[0]);
	int np = int(nnp[0]);
   	NumericVector prop0(1,0.0);
   	NumericVector prop1(1,0.0);
   	NumericVector datamin(1,0.0);
	NumericVector datamax(1,0.0);
   	NumericVector theta0old(1,0.0);
	NumericVector theta0new(1,0.0);
	NumericVector theta1old(1,0.0);
	NumericVector theta1new(1,0.0);
	NumericVector loglikdiff(1,0.0);
	NumericVector r(1,0.0);
	NumericVector uu(1,0.0);
	NumericVector vv(1,0.0);
	NumericVector postsamples0(M,0.0);
	NumericVector postsamples1(M,0.0);
	NumericVector logpost(M,0.0);
	NumericVector l0(1,0.0);
	NumericVector l1(1,0.0);
	NumericVector u0(1,0.0);
	NumericVector u1(1,0.0);
	theta0old(0) = bootmean0(0);
	theta1old(0) = bootmean1(0);
	NumericVector n0(1,0.0);
	NumericVector n1(1,0.0);
	NumericVector n2(1,0.0);
	NumericVector n0p(1,0.0);
	NumericVector n1p(1,0.0);
	NumericVector n2p(1,0.0);
	NumericVector F0_c0old(1,0.0);
	NumericVector F0_c0new(1,0.0);
	NumericVector F1_c0old(1,0.0);
	NumericVector F1_c0new(1,0.0);
	NumericVector F1_c1old(1,0.0);
	NumericVector F1_c1new(1,0.0);
	NumericVector F2_c1old(1,0.0);
	NumericVector F2_c1new(1,0.0);
	NumericVector F0_c0oldp(1,0.0);
	NumericVector F0_c0newp(1,0.0);
	NumericVector F1_c0oldp(1,0.0);
	NumericVector F1_c0newp(1,0.0);
	NumericVector F1_c1oldp(1,0.0);
	NumericVector F1_c1newp(1,0.0);
	NumericVector F2_c1oldp(1,0.0);
	NumericVector F2_c1newp(1,0.0);
	NumericVector YI(M,0.0);
	NumericVector YIl(1,0.0);
	NumericVector YIu(1,0.0);
	datamin(0) = data(0,1);
	datamax(0) = data(0,1);	
	NumericVector acc(1,0.0);
	
	for(int j=0; j<sN; j++){
		if(w(0)<=priorSched(j,0)){
			prop0(0) = priorSched(j,1);	
		}
		if(w(0)<=priorSched(j,2)){
			prop1(0) = priorSched(j,3);	
		}		
	}
	
	
	for(int k=0; k<n; k++){
		if(data(k,0)==1){
			n0(0) = n0(0) + 1.0;
			if(data(k,1)<=theta0old(0)){
				F0_c0old(0) = 	F0_c0old(0)+1.0;
			}	
		}
		else if(data(k,0)==2){
			n1(0) = n1(0) + 1.0;
			if(data(k,1)<=theta0old(0)){
				F1_c0old(0) = 	F1_c0old(0)+1.0;
			}
			if(data(k,1)<=theta1old(0)){
				F1_c1old(0) = 	F1_c1old(0)+1.0;
			}	
		}
		else {
			n2(0) = n2(0) + 1.0;
			if(data(k,1)<=theta1old(0)){
				F2_c1old(0) = 	F2_c1old(0)+1.0;
			}
		}
		if(datamin(0) > data(k,1)){
			datamin(0) = data(k,1);	
		}
		if(datamax(0) < data(k,1)){
			datamax(0) = data(k,1);	
		}
	}
	F0_c0old(0) = 	F0_c0old(0)/n0(0);
	F1_c0old(0) = 	F1_c0old(0)/n1(0);
	F1_c1old(0) = 	F1_c1old(0)/n1(0);
	F2_c1old(0) = 	F2_c1old(0)/n2(0);
	if(priorweight[0]>0.0){
	for(int k=0; k<np; k++){
		if(priordata(k,0)==1){
			n0p(0) = n0p(0) + 1.0;
			if(priordata(k,1)<=theta0old(0)){
				F0_c0oldp(0) = 	F0_c0oldp(0)+1.0;
			}	
		}
		else if(priordata(k,0)==2){
			n1p(0) = n1p(0) + 1.0;
			if(priordata(k,1)<=theta0old(0)){
				F1_c0oldp(0) = 	F1_c0oldp(0)+1.0;
			}
			if(priordata(k,1)<=theta1old(0)){
				F1_c1oldp(0) = 	F1_c1oldp(0)+1.0;
			}	
		}
		else {
			n2p(0) = n2p(0) + 1.0;
			if(priordata(k,1)<=theta1old(0)){
				F2_c1oldp(0) = 	F2_c1oldp(0)+1.0;
			}
		}
	}
	F0_c0oldp(0) = 	F0_c0oldp(0)/n0p(0);
	F1_c0oldp(0) = 	F1_c0oldp(0)/n1p(0);
	F1_c1oldp(0) = 	F1_c1oldp(0)/n1p(0);
	F2_c1oldp(0) = 	F2_c1oldp(0)/n2p(0);}
	
	for(int j=0; j<M; j++) {
		theta0new(0) = R::rnorm(theta0old(0), prop0(0));
		theta1new(0) = R::rnorm(theta1old(0), prop1(0));
		
		
		if( (theta0new(0)>theta1new(0)) ||  (theta0new(0)<datamin(0)) || (theta1new(0)>datamax(0))   ){
			postsamples0(j) = theta0old(0);
			postsamples1(j) = theta1old(0);
			logpost(j) = w[0]*n*(F0_c0old(0)+F1_c1old(0)-F1_c0old(0)-F2_c1old(0))+w[0]*priorweight[0]*np*(F0_c0oldp(0)+F1_c1oldp(0)-F1_c0oldp(0)-F2_c1oldp(0));
			YI(j) = F0_c0old(0) + F1_c1old(0) - F1_c0old(0) - F2_c1old(0);
		}else {
		
		loglikdiff(0) = 0.0;
		for(int k=0; k<n; k++){
			if(data(k,0)==1){
				if(data(k,1)<=theta0new(0)){
					F0_c0new(0) = 	F0_c0new(0)+1.0;
				}	
			}
			else if(data(k,0)==2){
				if(data(k,1)<=theta0new(0)){
					F1_c0new(0) = 	F1_c0new(0)+1.0;
				}
				if(data(k,1)<=theta1new(0)){
					F1_c1new(0) = 	F1_c1new(0)+1.0;
				}	
			}
			else {
				if(data(k,1)<=theta1new(0)){
					F2_c1new(0) = 	F2_c1new(0)+1.0;
				}
			}
		}
		F0_c0new(0) = 	F0_c0new(0)/n0(0);
		F1_c0new(0) = 	F1_c0new(0)/n1(0);
		F1_c1new(0) = 	F1_c1new(0)/n1(0);
		F2_c1new(0) = 	F2_c1new(0)/n2(0);
		if(priorweight[0]>0.0){
		for(int k=0; k<np; k++){
			if(priordata(k,0)==1){
				if(priordata(k,1)<=theta0new(0)){
					F0_c0newp(0) = 	F0_c0newp(0)+1.0;
				}	
			}
			else if(priordata(k,0)==2){
				if(priordata(k,1)<=theta0new(0)){
					F1_c0newp(0) = 	F1_c0newp(0)+1.0;
				}
				if(priordata(k,1)<=theta1new(0)){
					F1_c1newp(0) = 	F1_c1newp(0)+1.0;
				}	
			}
			else {
				if(priordata(k,1)<=theta1new(0)){
					F2_c1newp(0) = 	F2_c1newp(0)+1.0;
				}
			}
		}
		F0_c0newp(0) = 	F0_c0newp(0)/n0p(0);
		F1_c0newp(0) = 	F1_c0newp(0)/n1p(0);
		F1_c1newp(0) = 	F1_c1newp(0)/n1p(0);
		F2_c1newp(0) = 	F2_c1newp(0)/n2p(0);}
	
		loglikdiff(0) = w[0]*n*((F0_c0new(0)+F1_c1new(0)-F1_c0new(0)-F2_c1new(0))-(F0_c0old(0)+F1_c1old(0)-F1_c0old(0)-F2_c1old(0)))  +  w[0]*priorweight[0]*np*((F0_c0newp(0)+F1_c1newp(0)-F1_c0newp(0)-F2_c1newp(0))-(F0_c0oldp(0)+F1_c1oldp(0)-F1_c0oldp(0)-F2_c1oldp(0)));
		loglikdiff(0) = fmin(std::exp(loglikdiff(0)), 1.0);
		uu[0] = R::runif(0.0,1.0);
		if(uu(0) <= loglikdiff(0)) {
			postsamples0(j) = theta0new(0);
			postsamples1(j) = theta1new(0);
			logpost(j) = w[0]*n*(F0_c0new(0)+F1_c1new(0)-F1_c0new(0)-F2_c1new(0))+w[0]*priorweight[0]*np*(F0_c0newp(0)+F1_c1newp(0)-F1_c0newp(0)-F2_c1newp(0));
			theta0old(0) = theta0new(0);
			theta1old(0) = theta1new(0);
			F0_c0old(0)=F0_c0new(0);
			F1_c0old(0)=F1_c0new(0);
			F1_c1old(0)=F1_c1new(0);
			F2_c1old(0)=F2_c1new(0);
			F0_c0oldp(0)=F0_c0newp(0);
			F1_c0oldp(0)=F1_c0newp(0);
			F1_c1oldp(0)=F1_c1newp(0);
			F2_c1oldp(0)=F2_c1newp(0);
			YI(j) = F0_c0new(0) + F1_c1new(0) - F1_c0new(0) - F2_c1new(0);
			acc(0) = acc(0) + 1.0;
		}
		else {
			postsamples0(j) = theta0old(0);
			postsamples1(j) = theta1old(0);
			logpost(j) = w[0]*n*(F0_c0old(0)+F1_c1old(0)-F1_c0old(0)-F2_c1old(0))+w[0]*priorweight[0]*np*(F0_c0oldp(0)+F1_c1oldp(0)-F1_c0oldp(0)-F2_c1oldp(0));
			YI(j) = F0_c0old(0) + F1_c1old(0) - F1_c0old(0) - F2_c1old(0);
		}
		}
	}

	double templogpost;
	double tempYI;
	double temppost0;
	double temppost1;
	bool swapped;	
        for (int j = 0; j < M-1; j++){ 
		 swapped = false; 
       		 for (int k = 0; k < M-j-1; k++){  
           		 if (logpost(k) > logpost(k+1)){ 
				 swapped = true;
				 templogpost = logpost(k);
				 logpost(k) = logpost(k+1);
				 logpost(k+1) = templogpost;
				 tempYI = YI(k);
				 YI(k) = YI(k+1);
				 YI(k+1) = tempYI;
				 temppost0 = postsamples0(k);
				 postsamples0(k) = postsamples0(k+1);
				 postsamples0(k+1) = temppost0;				 
				 temppost1 = postsamples1(k);
				 postsamples1(k) = postsamples1(k+1);
				 postsamples1(k+1) = temppost1;
			 }
		 }
		 if (swapped == false) 
			break;
        }
			
	l0[0] = 100000;
	u0[0] = -100000;
	l1[0] = 100000;
	u1[0] = -100000;
	YIl[0] = 100000;
	YIu[0] = -100000;
	for(int i = (0.05*M-1); i<M; i++){
		if(l0[0]>postsamples0(i)){
			l0[0]=postsamples0(i);
		}
		if(u0[0]<postsamples0(i)){
			u0[0]=postsamples0(i);
		}
		if(l1[0]>postsamples1(i)){
			l1[0]=postsamples1(i);
		}
		if(u1[0]<postsamples1(i)){
			u1[0]=postsamples1(i);
		}
		if(YIl[0]>YI(i)){
			YIl[0]=YI(i);
		}
		if(YIu[0]<YI(i)){
			YIu[0]=YI(i);
		}
	}

	acc(0) = acc(0)/M;
	result = Rcpp::List::create(Rcpp::Named("l0") = l0,Rcpp::Named("u0") = u0,Rcpp::Named("l1") = l1,Rcpp::Named("u1") = u1,Rcpp::Named("YI") = YI,Rcpp::Named("YIl") = YIl,Rcpp::Named("YIu") = YIu,Rcpp::Named("datamax") = datamax,Rcpp::Named("datamin") = datamin, Rcpp::Named("acceptance_rate") = acc, Rcpp::Named("samples0") = postsamples0, Rcpp::Named("samples1") = postsamples1, Rcpp::Named("logpost") = logpost);

	return result;
}


Rcpp::List GibbsMCMCkde2(NumericVector nn, NumericMatrix data, NumericVector nnp, NumericMatrix priordata, NumericVector priorweight, NumericMatrix thetaboot,
	NumericVector bootmean0, NumericVector bootmean1, NumericVector kdecdflen, NumericMatrix kdecdfboot1, NumericMatrix kdecdfboot2, NumericMatrix kdecdfboot3, NumericMatrix kdecdfboot1p, NumericMatrix kdecdfboot2p, NumericMatrix kdecdfboot3p, NumericVector scheduleLen, NumericMatrix priorSched, NumericVector alpha, NumericVector M_samp, NumericVector w) {
   	
	List result;
	int M = int(M_samp[0]);
	int n = int(nn[0]);
	int sN = int(scheduleLen[0]);
	int np = int(nnp[0]);
	int kdeN = int(kdecdflen[0]);
   	NumericVector prop0(1,0.0);
   	NumericVector prop1(1,0.0);
   	NumericVector theta0old(1,0.0);
	NumericVector theta0new(1,0.0);
	NumericVector theta1old(1,0.0);
	NumericVector theta1new(1,0.0);
	NumericVector loglikdiff(1,0.0);
	NumericVector r(1,0.0);
	NumericVector uu(1,0.0);
	NumericVector vv(1,0.0);
	NumericVector postsamples0(M,0.0);
	NumericVector postsamples1(M,0.0);
	NumericVector logpost(M,0.0);
	NumericVector l0(1,0.0);
	NumericVector l1(1,0.0);
	NumericVector u0(1,0.0);
	NumericVector u1(1,0.0);
	theta0old(0) = bootmean0(0);
	theta1old(0) = bootmean1(0);
	NumericVector F0_c0old(1,0.0);
	NumericVector F0_c0new(1,0.0);
	NumericVector F1_c0old(1,0.0);
	NumericVector F1_c0new(1,0.0);
	NumericVector F1_c1old(1,0.0);
	NumericVector F1_c1new(1,0.0);
	NumericVector F2_c1old(1,0.0);
	NumericVector F2_c1new(1,0.0);
	NumericVector F0_c0oldp(1,0.0);
	NumericVector F0_c0newp(1,0.0);
	NumericVector F1_c0oldp(1,0.0);
	NumericVector F1_c0newp(1,0.0);
	NumericVector F1_c1oldp(1,0.0);
	NumericVector F1_c1newp(1,0.0);
	NumericVector F2_c1oldp(1,0.0);
	NumericVector F2_c1newp(1,0.0);
	NumericVector YI(M,0.0);
	NumericVector YIl(1,0.0);
	NumericVector YIu(1,0.0);
	NumericVector acc(1,0.0);
	
	for(int j=0; j<sN; j++){
		if(w(0)<=priorSched(j,0)){
			prop0(0) = priorSched(j,1);	
		}
		if(w(0)<=priorSched(j,2)){
			prop1(0) = priorSched(j,3);	
		}		
	}
	
	
	
	int k = 0;
	while(k < kdeN){
		if(theta0old(0)<=kdecdfboot1(k,0)){
			F0_c0old(0) = kdecdfboot1(k,1);
			break;
		}	
		k = k + 1;
	}
	k = 0;
	while(k < kdeN){
		if(theta0old(0)<=kdecdfboot2(k,0)){
			F1_c0old(0) = kdecdfboot2(k,1);
			break;
		}
		k = k + 1;
	}
	k = 0;
	while(k < kdeN){
		if(theta1old(0)<=kdecdfboot2(k,0)){
			F1_c1old(0) = kdecdfboot2(k,1);
			break;
		}
		k = k + 1;
	}
	k = 0;
	while(k < kdeN){
		if(theta1old(0)<=kdecdfboot3(k,0)){
			F2_c1old(0) = kdecdfboot3(k,1);
			break;
		}
		k = k + 1;
	}	
	
	
	if(priorweight(0)>0.0){
		k = 0;
		while(k < kdeN){
			if(theta0old(0)<=kdecdfboot1p(k,0)){
				F0_c0oldp(0) = kdecdfboot1p(k,1);
				break;
			}
			k = k + 1;
		}
		k = 0;
		while(k < kdeN){
			if(theta0old(0)<=kdecdfboot2p(k,0)){
				F1_c0oldp(0) = kdecdfboot2p(k,1);
				break;
			}	
			k = k + 1;
		}
		k = 0;
		while(k < kdeN){
			if(theta1old(0)<=kdecdfboot2p(k,0)){
				F1_c1oldp(0) = kdecdfboot2p(k,1);
				break;
			}
			k = k + 1;
		}
		k = 0;
		while(k < kdeN){
			if(theta1old(0)<=kdecdfboot3p(k,0)){
				F2_c1oldp(0) = kdecdfboot3p(k,1);
				break;
			}
			k = k + 1;
		}	
	}
	
	
	for(int j=0; j<M; j++) {
		theta0new(0) = R::rnorm(theta0old(0), prop0(0));
		theta1new(0) = R::rnorm(theta1old(0), prop1(0));
		
		
		if( (theta0new(0)>theta1new(0)) || (theta0new(0) < kdecdfboot1(0,0)) || (theta1new(0) > kdecdfboot3(kdeN-1,0))   ){
			postsamples0(j) = theta0old(0);
			postsamples1(j) = theta1old(0);
			logpost(j) = w[0]*n*(F0_c0old(0)+F1_c1old(0)-F1_c0old(0)-F2_c1old(0))+w[0]*priorweight[0]*np*(F0_c0oldp(0)+F1_c1oldp(0)-F1_c0oldp(0)-F2_c1oldp(0));
			YI(j) = F0_c0old(0) + F1_c1old(0) - F1_c0old(0) - F2_c1old(0);
		}else {
		
		loglikdiff(0) = 0.0;

		k = 0;
		while(k < kdeN){
			if(theta0new(0)<=kdecdfboot1(k,0)){
				F0_c0new(0) = kdecdfboot1(k,1);
				break;
			}
			k = k + 1;
		}
		k = 0;
		while(k < kdeN){
			if(theta0new(0)<=kdecdfboot2(k,0)){
				F1_c0new(0) = kdecdfboot2(k,1);
				break;
			}
			k = k + 1;
		}
		k = 0;
		while(k < kdeN){
			if(theta1new(0)<=kdecdfboot2(k,0)){
				F1_c1new(0) = kdecdfboot2(k,1);
				break;
			}
			k = k + 1;
		}
		k = 0;
		while(k < kdeN){
			if(theta1new(0)<=kdecdfboot3(k,0)){
				F2_c1new(0) = kdecdfboot3(k,1);
				break;
			}
			k = k + 1;
		}	
	
	
		if(priorweight(0)>0.0){
			k = 0;
			while(k < kdeN){
				if(theta0new(0)<=kdecdfboot1p(k,0)){
					F0_c0newp(0) = kdecdfboot1p(k,1);
					break;
				}
				k = k + 1;
			}
			k = 0;
			while(k < kdeN){
				if(theta0new(0)<=kdecdfboot2p(k,0)){
					F1_c0newp(0) = kdecdfboot2p(k,1);
					break;
				}
				k = k + 1;
			}
			k = 0;
			while(k < kdeN){
				if(theta1new(0)<=kdecdfboot2p(k,0)){
					F1_c1newp(0) = kdecdfboot2p(k,1);
					break;
				}	
				k = k + 1;
			}
			k = 0;
			while(k < kdeN){
				if(theta1new(0)<=kdecdfboot3p(k,0)){
					F2_c1newp(0) = kdecdfboot3p(k,1);
					break;
				}	
				k = k + 1;
			}	
		}
	
		loglikdiff(0) = w[0]*n*((F0_c0new(0)+F1_c1new(0)-F1_c0new(0)-F2_c1new(0))-(F0_c0old(0)+F1_c1old(0)-F1_c0old(0)-F2_c1old(0)))  +  w[0]*priorweight[0]*np*((F0_c0newp(0)+F1_c1newp(0)-F1_c0newp(0)-F2_c1newp(0))-(F0_c0oldp(0)+F1_c1oldp(0)-F1_c0oldp(0)-F2_c1oldp(0)));
		loglikdiff(0) = fmin(std::exp(loglikdiff(0)), 1.0);
		uu[0] = R::runif(0.0,1.0);
		if(uu(0) <= loglikdiff(0)) {
			postsamples0(j) = theta0new(0);
			postsamples1(j) = theta1new(0);
			logpost(j) = w[0]*n*(F0_c0new(0)+F1_c1new(0)-F1_c0new(0)-F2_c1new(0))+w[0]*priorweight[0]*np*(F0_c0newp(0)+F1_c1newp(0)-F1_c0newp(0)-F2_c1newp(0));
			theta0old(0) = theta0new(0);
			theta1old(0) = theta1new(0);
			F0_c0old(0)=F0_c0new(0);
			F1_c0old(0)=F1_c0new(0);
			F1_c1old(0)=F1_c1new(0);
			F2_c1old(0)=F2_c1new(0);
			F0_c0oldp(0)=F0_c0newp(0);
			F1_c0oldp(0)=F1_c0newp(0);
			F1_c1oldp(0)=F1_c1newp(0);
			F2_c1oldp(0)=F2_c1newp(0);
			YI(j) = F0_c0new(0) + F1_c1new(0) - F1_c0new(0) - F2_c1new(0);
			acc(0) = acc(0) + 1.0;
		}
		else {
			postsamples0(j) = theta0old(0);
			postsamples1(j) = theta1old(0);
			logpost(j) = w[0]*n*(F0_c0old(0)+F1_c1old(0)-F1_c0old(0)-F2_c1old(0))+w[0]*priorweight[0]*np*(F0_c0oldp(0)+F1_c1oldp(0)-F1_c0oldp(0)-F2_c1oldp(0));
			YI(j) = F0_c0old(0) + F1_c1old(0) - F1_c0old(0) - F2_c1old(0);
		}
		}
	}

	double templogpost;
	double tempYI;
	double temppost0;
	double temppost1;
	bool swapped;	
        for (int j = 0; j < M-1; j++){ 
		 swapped = false; 
       		 for (int k = 0; k < M-j-1; k++){  
           		 if (logpost(k) > logpost(k+1)){ 
				 swapped = true;
				 templogpost = logpost(k);
				 logpost(k) = logpost(k+1);
				 logpost(k+1) = templogpost;
				 tempYI = YI(k);
				 YI(k) = YI(k+1);
				 YI(k+1) = tempYI;
				 temppost0 = postsamples0(k);
				 postsamples0(k) = postsamples0(k+1);
				 postsamples0(k+1) = temppost0;				 
				 temppost1 = postsamples1(k);
				 postsamples1(k) = postsamples1(k+1);
				 postsamples1(k+1) = temppost1;
			 }
		 }
		 if (swapped == false) 
			break;
        }
			
	l0[0] = 100000;
	u0[0] = -100000;
	l1[0] = 100000;
	u1[0] = -100000;
	YIl[0] = 100000;
	YIu[0] = -100000;
	for(int i = (0.05*M-1); i<M; i++){
		if(l0[0]>postsamples0(i)){
			l0[0]=postsamples0(i);
		}
		if(u0[0]<postsamples0(i)){
			u0[0]=postsamples0(i);
		}
		if(l1[0]>postsamples1(i)){
			l1[0]=postsamples1(i);
		}
		if(u1[0]<postsamples1(i)){
			u1[0]=postsamples1(i);
		}
		if(YIl[0]>YI(i)){
			YIl[0]=YI(i);
		}
		if(YIu[0]<YI(i)){
			YIu[0]=YI(i);
		}
	}

	acc(0) = acc(0)/M;
	result = Rcpp::List::create(Rcpp::Named("l0") = l0,Rcpp::Named("u0") = u0,Rcpp::Named("l1") = l1,Rcpp::Named("u1") = u1,Rcpp::Named("YI") = YI,Rcpp::Named("YIl") = YIl,Rcpp::Named("YIu") = YIu, Rcpp::Named("acceptance_rate") = acc, Rcpp::Named("samples0") = postsamples0, Rcpp::Named("samples1") = postsamples1, Rcpp::Named("logpost") = logpost);

	return result;
}



/* adaptive proposal
// [[Rcpp::export]]
Rcpp::List GibbsMCMC2(NumericVector nn, NumericMatrix data, NumericVector nnp, NumericMatrix priordata, NumericVector priorweight, NumericMatrix thetaboot,
	NumericVector bootmean0, NumericVector bootmean1, NumericVector alpha, NumericVector M_samp, NumericVector w) {
   	
	List result;
	int M = int(M_samp[0]);
	int n = int(nn[0]);
	int np = int(nnp[0]);
   	NumericVector datamin(1,0.0);
	NumericVector datamax(1,0.0);
   	NumericVector theta0old(1,0.0);
	NumericVector theta0new(1,0.0);
	NumericVector theta1old(1,0.0);
	NumericVector theta1new(1,0.0);
	NumericVector loglikdiff(1,0.0);
	NumericVector r(1,0.0);
	NumericVector uu(1,0.0);
	NumericVector vv(1,0.0);
	NumericVector postsamples0(M,0.0);
	NumericVector postsamples1(M,0.0);
	NumericVector l0(1,0.0);
	NumericVector l1(1,0.0);
	NumericVector u0(1,0.0);
	NumericVector u1(1,0.0);
	theta0old(0) = bootmean0(0);
	theta1old(0) = bootmean1(0);
	NumericVector n0(1,0.0);
	NumericVector n1(1,0.0);
	NumericVector n2(1,0.0);
	NumericVector n0p(1,0.0);
	NumericVector n1p(1,0.0);
	NumericVector n2p(1,0.0);
	NumericVector F0_c0old(1,0.0);
	NumericVector F0_c0new(1,0.0);
	NumericVector F1_c0old(1,0.0);
	NumericVector F1_c0new(1,0.0);
	NumericVector F1_c1old(1,0.0);
	NumericVector F1_c1new(1,0.0);
	NumericVector F2_c1old(1,0.0);
	NumericVector F2_c1new(1,0.0);
	NumericVector F0_c0oldp(1,0.0);
	NumericVector F0_c0newp(1,0.0);
	NumericVector F1_c0oldp(1,0.0);
	NumericVector F1_c0newp(1,0.0);
	NumericVector F1_c1oldp(1,0.0);
	NumericVector F1_c1newp(1,0.0);
	NumericVector F2_c1oldp(1,0.0);
	NumericVector F2_c1newp(1,0.0);
	NumericVector sumsamp0(1,0.0);
	NumericVector sumsamp1(1,0.0);
	NumericVector sumsamp0sq(1,0.0);
	NumericVector sumsamp1sq(1,0.0);
	NumericVector sumsamp01(1,0.0);
	NumericVector s2x(1,0.0);
	NumericVector s2y(1,0.0);
	NumericVector sxy(1,0.0);
	NumericVector YI(M,0.0);
	NumericVector YIl(1,0.0);
	NumericVector YIu(1,0.0);
	datamin(0) = data(0,1);
	datamax(0) = data(0,1);	
	NumericVector acc(1,0.0);
	
	
	for(int k=0; k<n; k++){
		if(data(k,0)==1){
			n0(0) = n0(0) + 1.0;
			if(data(k,1)<=theta0old(0)){
				F0_c0old(0) = 	F0_c0old(0)+1.0;
			}	
		}
		else if(data(k,0)==2){
			n1(0) = n1(0) + 1.0;
			if(data(k,1)<=theta0old(0)){
				F1_c0old(0) = 	F1_c0old(0)+1.0;
			}
			if(data(k,1)<=theta1old(0)){
				F1_c1old(0) = 	F1_c1old(0)+1.0;
			}	
		}
		else {
			n2(0) = n2(0) + 1.0;
			if(data(k,1)<=theta1old(0)){
				F2_c1old(0) = 	F2_c1old(0)+1.0;
			}
		}
		if(datamin(0) > data(k,1)){
			datamin(0) = data(k,1);	
		}
		if(datamax(0) < data(k,1)){
			datamax(0) = data(k,1);	
		}
	}
	F0_c0old(0) = 	F0_c0old(0)/n0(0);
	F1_c0old(0) = 	F1_c0old(0)/n1(0);
	F1_c1old(0) = 	F1_c1old(0)/n1(0);
	F2_c1old(0) = 	F2_c1old(0)/n2(0);
	
	for(int k=0; k<np; k++){
		if(priordata(k,0)==1){
			n0p(0) = n0p(0) + 1.0;
			if(priordata(k,1)<=theta0old(0)){
				F0_c0oldp(0) = 	F0_c0oldp(0)+1.0;
			}	
		}
		else if(priordata(k,0)==2){
			n1p(0) = n1p(0) + 1.0;
			if(priordata(k,1)<=theta0old(0)){
				F1_c0oldp(0) = 	F1_c0oldp(0)+1.0;
			}
			if(priordata(k,1)<=theta1old(0)){
				F1_c1oldp(0) = 	F1_c1oldp(0)+1.0;
			}	
		}
		else {
			n2p(0) = n2p(0) + 1.0;
			if(priordata(k,1)<=theta1old(0)){
				F2_c1oldp(0) = 	F2_c1oldp(0)+1.0;
			}
		}
	}
	F0_c0oldp(0) = 	F0_c0oldp(0)/n0p(0);
	F1_c0oldp(0) = 	F1_c0oldp(0)/n1p(0);
	F1_c1oldp(0) = 	F1_c1oldp(0)/n1p(0);
	F2_c1oldp(0) = 	F2_c1oldp(0)/n2p(0);
	
	for(int j=0; j<M; j++) {
		if(j<=16){
			theta0new(0) = R::rnorm(theta0old(0), .1);
			theta1new(0) = R::rnorm(theta1old(0), .1);
		}
		else {
			vv[0] = R::runif(0.0,1.0);
			if(vv[0]<=0.95){
				theta0new(0) = R::rnorm(0.0, 1.0);
				theta1new(0) = R::rnorm(0.0, 1.0);
				s2x(0) = (2.38*(2.38/2.0))*((sumsamp0sq(0)/j) - pow((sumsamp0(0)/j),2.0));
				s2y(0) = (2.38*(2.38/2.0))*((sumsamp1sq(0)/j) - pow((sumsamp1(0)/j),2.0));
				sxy(0) = (2.38*(2.38/2.0))*((sumsamp01(0)/j) - (sumsamp0(0)*sumsamp1(0)/(j*j)));
				theta0new(0) = (theta0new(0)*sqrt(s2x(0)))+theta0old(0);
				theta1new(0) = (theta1new(0)*sqrt(s2y(0)-(pow(sxy(0),2.0)*(1/s2x(0))))) + (theta0new(0)*sxy(0)*(1/sqrt(s2x(0))))+theta1old(0);
			}
			else {
				theta0new(0) = R::rnorm(theta0old(0), 0.1);
				theta1new(0) = R::rnorm(theta1old(0), 0.1);					
			}
		}
		if(theta0new(0)>theta1new(0)){
			theta0new(0) = theta0old(0);
			theta1new(0) = theta1old(0);			
		}
		if(theta0new(0)<datamin(0)){
			theta0new(0) = theta0old(0);			
		}
		if(theta1new(0)>datamax(0)){
			theta1new(0) = theta1old(0);			
		}
		loglikdiff(0) = 0.0;
		for(int k=0; k<n; k++){
			if(data(k,0)==1){
				if(data(k,1)<=theta0new(0)){
					F0_c0new(0) = 	F0_c0new(0)+1.0;
				}	
			}
			else if(data(k,0)==2){
				if(data(k,1)<=theta0new(0)){
					F1_c0new(0) = 	F1_c0new(0)+1.0;
				}
				if(data(k,1)<=theta1new(0)){
					F1_c1new(0) = 	F1_c1new(0)+1.0;
				}	
			}
			else {
				if(data(k,1)<=theta1new(0)){
					F2_c1new(0) = 	F2_c1new(0)+1.0;
				}
			}
		}
		F0_c0new(0) = 	F0_c0new(0)/n0(0);
		F1_c0new(0) = 	F1_c0new(0)/n1(0);
		F1_c1new(0) = 	F1_c1new(0)/n1(0);
		F2_c1new(0) = 	F2_c1new(0)/n2(0);
		
		for(int k=0; k<np; k++){
			if(priordata(k,0)==1){
				if(priordata(k,1)<=theta0new(0)){
					F0_c0newp(0) = 	F0_c0newp(0)+1.0;
				}	
			}
			else if(priordata(k,0)==2){
				if(priordata(k,1)<=theta0new(0)){
					F1_c0newp(0) = 	F1_c0newp(0)+1.0;
				}
				if(priordata(k,1)<=theta1new(0)){
					F1_c1newp(0) = 	F1_c1newp(0)+1.0;
				}	
			}
			else {
				if(priordata(k,1)<=theta1new(0)){
					F2_c1newp(0) = 	F2_c1newp(0)+1.0;
				}
			}
		}
		F0_c0newp(0) = 	F0_c0newp(0)/n0p(0);
		F1_c0newp(0) = 	F1_c0newp(0)/n1p(0);
		F1_c1newp(0) = 	F1_c1newp(0)/n1p(0);
		F2_c1newp(0) = 	F2_c1newp(0)/n2p(0);
	
		loglikdiff(0) = w[0]*n*((F0_c0new(0)+F1_c1new(0)-F1_c0new(0)-F2_c1new(0))-(F0_c0old(0)+F1_c1old(0)-F1_c0old(0)-F2_c1old(0)))  +  w[0]*priorweight[0]*np*((F0_c0newp(0)+F1_c1newp(0)-F1_c0newp(0)-F2_c1newp(0))-(F0_c0oldp(0)+F1_c1oldp(0)-F1_c0oldp(0)-F2_c1oldp(0)));
		loglikdiff(0) = fmin(std::exp(loglikdiff(0)), 1.0);
		uu[0] = R::runif(0.0,1.0);
		if(uu(0) <= loglikdiff(0)) {
			postsamples0(j) = theta0new(0);
			postsamples1(j) = theta1new(0);
			theta0old(0) = theta0new(0);
			theta1old(0) = theta1new(0);
			F0_c0old(0)=F0_c0new(0);
			F1_c0old(0)=F1_c0new(0);
			F1_c1old(0)=F1_c1new(0);
			F2_c1old(0)=F2_c1new(0);
			F0_c0oldp(0)=F0_c0newp(0);
			F1_c0oldp(0)=F1_c0newp(0);
			F1_c1oldp(0)=F1_c1newp(0);
			F2_c1oldp(0)=F2_c1newp(0);
			sumsamp0(0)=sumsamp0(0)+theta0new(0);
			sumsamp0sq(0)=sumsamp0sq(0)+theta0new(0)*theta0new(0);
			sumsamp1(0)=sumsamp1(0)+theta1new(0);
			sumsamp1sq(0)=sumsamp1sq(0)+theta1new(0)*theta1new(0);
			sumsamp01(0)=sumsamp01(0)+theta1new(0)*theta0new(0);
			YI(j) = F0_c0new(0) + F1_c1new(0) - F1_c0new(0) - F2_c1new(0);
			acc(0) = acc(0) + 1.0;
		}
		else {
			postsamples0(j) = theta0old(0);
			postsamples1(j) = theta1old(0);
			sumsamp0(0)=sumsamp0(0)+theta0old(0);
			sumsamp0sq(0)=sumsamp0sq(0)+theta0old(0)*theta0old(0);
			sumsamp1(0)=sumsamp1(0)+theta1old(0);
			sumsamp1sq(0)=sumsamp1sq(0)+theta1old(0)*theta1old(0);
			sumsamp01(0)=sumsamp01(0)+theta1old(0)*theta0old(0);
			YI(j) = F0_c0old(0) + F1_c1old(0) - F1_c0old(0) - F2_c1old(0);
		}
	}
	
	double tempYI;
	double temppost0;
	double temppost1;
	bool swapped;	
        for (int i = 0; i < M-1; i++){ 
		 swapped = false; 
       		 for (int j = 0; j < M-i-1; j++){  
           		 if (YI(j) > YI(j+1)){ 
				 swapped = true;
				 tempYI = YI(j);
				 YI(j) = YI(j+1);
				 YI(j+1) = tempYI;
				 temppost0 = postsamples0(j);
				 postsamples0(j) = postsamples0(j+1);
				 postsamples0(j+1) = temppost0;				 
				 temppost1 = postsamples1(j);
				 postsamples1(j) = postsamples1(j+1);
				 postsamples1(j+1) = temppost1;
			 }
		 }
		 if (swapped == false) 
			break;
        }
			
	l0[0] = 100000;
	u0[0] = -100000;
	l1[0] = 100000;
	u1[0] = -100000;
	for(int i = (0.05*M-1); i<M; i++){
		if(l0[0]>postsamples0(i)){
			l0[0]=postsamples0(i);
		}
		if(u0[0]<postsamples0(i)){
			u0[0]=postsamples0(i);
		}
		if(l1[0]>postsamples1(i)){
			l1[0]=postsamples1(i);
		}
		if(u1[0]<postsamples1(i)){
			u1[0]=postsamples1(i);
		}
	}
	YIl[0] = YI(0.025*M-1);
	YIu[0] = YI(0.975*M-1);
	acc(0) = acc(0)/M;
	result = Rcpp::List::create(Rcpp::Named("l0") = l0,Rcpp::Named("u0") = u0,Rcpp::Named("l1") = l1,Rcpp::Named("u1") = u1,Rcpp::Named("YI") = YI,Rcpp::Named("YIl") = YIl,Rcpp::Named("YIu") = YIu,Rcpp::Named("datamax") = datamax,Rcpp::Named("datamin") = datamin, Rcpp::Named("acceptance_rate") = acc, Rcpp::Named("samples0") = postsamples0, Rcpp::Named("samples1") = postsamples1);

	return result;
}
*/
/* non-adaptive version
// [[Rcpp::export]]
Rcpp::List GibbsMCMC2(NumericVector nn, NumericMatrix data, NumericMatrix thetaboot,
	NumericVector bootmean0, NumericVector bootmean1, NumericVector alpha, NumericVector M_samp, NumericVector w) {
   	
	List result;
	int M = int(M_samp[0]);
	int n = int(nn[0]);
   	NumericVector theta0old(1,0.0);
	NumericVector theta0new(1,0.0);
	NumericVector theta1old(1,0.0);
	NumericVector theta1new(1,0.0);
	NumericVector loglikdiff(1,0.0);
	NumericVector r(1,0.0);
	NumericVector uu(1,0.0);
	NumericVector postsamples0(M,0.0);
	NumericVector postsamples1(M,0.0);
	NumericVector l0(1,0.0);
	NumericVector l1(1,0.0);
	NumericVector u0(1,0.0);
	NumericVector u1(1,0.0);
	NumericVector n0(1,0.0);
	NumericVector n1(1,0.0);
	NumericVector n2(1,0.0);
	NumericVector F0_c0old(1,0.0);
	NumericVector F0_c0new(1,0.0);
	NumericVector F1_c0old(1,0.0);
	NumericVector F1_c0new(1,0.0);
	NumericVector F1_c1old(1,0.0);
	NumericVector F1_c1new(1,0.0);
	NumericVector F2_c1old(1,0.0);
	NumericVector F2_c1new(1,0.0);
	theta0old = bootmean0;
	theta1old = bootmean1;
	
	for(int k=0; k<n; k++){
		if(data(k,0)==1){
			n0(0) = n0(0) + 1;
			if(data(k,1)<=theta0old(0)){
				F0_c0old(0) = 	F0_c0old(0)+1;
			}	
		}
		else if(data(k,0)==2){
			n1(0) = n1(0) + 1;
			if(data(k,1)<=theta0old(0)){
				F1_c0old(0) = 	F1_c0old(0)+1;
			}
			if(data(k,1)<=theta1old(0)){
				F1_c1old(0) = 	F1_c1old(0)+1;
			}	
		}
		else {
			n2(0) = n2(0) + 1;
			if(data(k,1)<=theta1old(0)){
				F2_c1old(0) = 	F2_c1old(0)+1;
			}
		}
	}
	F0_c0old(0) = 	F0_c0old(0)/n0(0);
	F1_c0old(0) = 	F1_c0old(0)/n1(0);
	F1_c1old(0) = 	F1_c1old(0)/n1(0);
	F2_c1old(0) = 	F2_c1old(0)/n2(0);
	
	for(int j=0; j<(M+100); j++) {
		theta0new(0) = R::rnorm(theta0old(0), 0.5);
		loglikdiff(0) = 0.0;
		for(int k=0; k<n; k++){
			if(data(k,0)==1){
				if(data(k,1)<=theta0new(0)){
					F0_c0new(0) = 	F0_c0new(0)+1;
				}	
			}
			else if(data(k,0)==2){
				if(data(k,1)<=theta0new(0)){
					F1_c0new(0) = 	F1_c0new(0)+1;
				}
			}
		}
		F0_c0new(0) = 	F0_c0new(0)/n0(0);
		F1_c0new(0) = 	F1_c0new(0)/n1(0);
		loglikdiff(0) = w[0]*10.0*(F0_c0new(0)-F0_c0old(0)-F1_c0new(0)+F1_c0old(0));
		loglikdiff(0) = fmin(std::exp(loglikdiff(0)), 1.0);
		uu[0] = R::runif(0.0,1.0);
      		if((uu(0) <= loglikdiff(0)) && (j>99)) {
			postsamples0(j-100) = theta0new(0);
			theta0old(0) = theta0new(0); 
			F0_c0old(0)=F0_c0new(0);
			F1_c0old(0)=F1_c0new(0);
      		}
		else if(j>99){
			postsamples0(j-100) = theta0old(0);	
		}
		theta1new[0] = R::rnorm(theta1old(0), 0.5);
		loglikdiff(0) = 0.0;
		for(int k=0; k<n; k++){
			if(data(k,0)==2){
				if(data(k,1)<=theta1new(0)){
					F1_c1new(0) = 	F1_c1new(0)+1;
				}	
			}
			else {
				if(data(k,1)<=theta1new(0)){
					F2_c1new(0) = 	F2_c1new(0)+1;
				}
			}
		}
		F1_c1new(0) = 	F1_c1new(0)/n1(0);
		F2_c1new(0) = 	F2_c1new(0)/n2(0);
		loglikdiff(0) = w[0]*10.0*(F1_c1new(0)-F1_c1old(0)-F2_c1new(0)+F2_c1old(0));		
		loglikdiff(0) = fmin(std::exp(loglikdiff(0)), 1.0);
		uu[0] = R::runif(0.0,1.0);
      		if((uu(0) <= loglikdiff(0)) && (j>99)) {
			postsamples1(j-100) = theta1new(0);
			theta1old(0) = theta1new(0); 
			F1_c1old(0)=F1_c1new(0);
			F2_c1old(0)=F2_c1new(0);
      		}
		else if(j>99){
			postsamples1(j-100) = theta1old(0);	
		}
	}
	std::sort(postsamples0.begin(), postsamples0.end());
	std::sort(postsamples1.begin(), postsamples1.end());
	l0[0] = postsamples0(0.025*M);
	u0[0] = postsamples0(0.975*M);
	l1[0] = postsamples1(0.025*M);
	u1[0] = postsamples1(0.975*M);
	
	result = Rcpp::List::create(Rcpp::Named("l0") = l0[0],Rcpp::Named("u0") = u0[0],Rcpp::Named("l1") = l1[0],Rcpp::Named("u1") = u1[0]);

	return result;
}

*/

struct GPCYI_yi_mcmc_parallel : public Worker {

	const RVector<double> nn;
	const RMatrix<double> data;
	const RVector<double> nnp;
	const RMatrix<double> priordata;
	const RVector<double> priorweight;
	const RMatrix<double> thetaboot;
	const RVector<double> bootmean0;
	const RVector<double> bootmean1;
	const RMatrix<double> databoot;
	const RVector<double> scheduleLen;
	const RMatrix<double> priorSched;
	const RVector<double> alpha;
	const RVector<double> M_samp;
	const RVector<double> B_resamp;
	const RVector<double> w;
	RVector<double> cover;

   // initialize with source and destination
   GPCYI_yi_mcmc_parallel(const NumericVector nn, const NumericMatrix data, const NumericVector nnp, const NumericMatrix priordata, const NumericVector priorweight, const NumericMatrix thetaboot,
	const NumericVector bootmean0, const NumericVector bootmean1, const NumericMatrix databoot, const NumericVector scheduleLen, const NumericMatrix priorSched,
	const NumericVector alpha, const NumericVector M_samp, const NumericVector B_resamp,
	const NumericVector w, NumericVector cover) 
			: nn(nn), data(data), nnp(nnp), priordata(priordata), priorweight(priorweight),  thetaboot(thetaboot), bootmean0(bootmean0), bootmean1(bootmean1), databoot(databoot), scheduleLen(scheduleLen), priorSched(priorSched), alpha(alpha), M_samp(M_samp), B_resamp(B_resamp), w(w), cover(cover) {}   

   // operator
void operator()(std::size_t begin, std::size_t end) {
		for (std::size_t i = begin; i < end; i++) {
			cover[i] = GibbsMCMC(nn, data, nnp, priordata, priorweight, thetaboot, bootmean0, bootmean1, databoot, scheduleLen, priorSched, alpha, M_samp, w, i);	
		}
	}
};

// [[Rcpp::export]]
NumericVector rcpp_parallel_yi(NumericVector nn, NumericMatrix data, NumericVector nnp, NumericMatrix priordata, NumericVector priorweight, NumericMatrix thetaboot, NumericVector bootmean0,
	NumericVector bootmean1, NumericMatrix databoot, NumericVector scheduleLen, NumericMatrix priorSched, NumericVector alpha, NumericVector M_samp, NumericVector B_resamp,
	NumericVector w) {
	
   int B = int(B_resamp[0]);
   // allocate the matrix we will return
   NumericVector cover(B,2.0); 

   // create the worker
   GPCYI_yi_mcmc_parallel gpcWorker(nn, data, nnp, priordata, priorweight, thetaboot, bootmean0, bootmean1, databoot, scheduleLen, priorSched, alpha, M_samp, B_resamp, w, cover);
     
   // call it with parallelFor
   
   parallelFor(0, B, gpcWorker);

   return cover;
}


struct GPCYI_yi_kde_mcmc_parallel : public Worker {

	const RVector<double> nn;
	const RMatrix<double> data;
	const RVector<double> nnp;
	const RMatrix<double> priordata;
	const RVector<double> priorweight;
	const RMatrix<double> thetaboot;
	const RVector<double> bootmean0;
	const RVector<double> bootmean1;
	const RVector<double> bootmeanYI;
	const RVector<double> kdecdflen;
	const RMatrix<double> kdecdfboot1;
	const RMatrix<double> kdecdfboot2;
	const RMatrix<double> kdecdfboot3;
	const RMatrix<double> kdecdfboot1p;
	const RMatrix<double> kdecdfboot2p;
	const RMatrix<double> kdecdfboot3p;
	const RVector<double> scheduleLen;
	const RMatrix<double> priorSched;
	const RVector<double> alpha;
	const RVector<double> M_samp;
	const RVector<double> B_resamp;
	const RVector<double> w;
	RVector<double> cover;

   // initialize with source and destination
   GPCYI_yi_kde_mcmc_parallel(const NumericVector nn, const NumericMatrix data, const NumericVector nnp, const NumericMatrix priordata, const NumericVector priorweight, const NumericMatrix thetaboot,
	const NumericVector bootmean0, const NumericVector bootmean1, const NumericVector bootmeanYI, const NumericVector kdecdflen, const NumericMatrix kdecdfboot1, const NumericMatrix kdecdfboot2, const NumericMatrix kdecdfboot3, const NumericMatrix kdecdfboot1p, const NumericMatrix kdecdfboot2p, const NumericMatrix kdecdfboot3p, const NumericVector scheduleLen, const NumericMatrix priorSched,
	const NumericVector alpha, const NumericVector M_samp, const NumericVector B_resamp,
	const NumericVector w, NumericVector cover) 
			: nn(nn), data(data), nnp(nnp), priordata(priordata), priorweight(priorweight),  thetaboot(thetaboot), bootmean0(bootmean0), bootmean1(bootmean1), bootmeanYI(bootmeanYI), kdecdflen(kdecdflen), kdecdfboot1(kdecdfboot1), kdecdfboot2(kdecdfboot2), kdecdfboot3(kdecdfboot3), kdecdfboot1p(kdecdfboot1p), kdecdfboot2p(kdecdfboot2p), kdecdfboot3p(kdecdfboot3p), scheduleLen(scheduleLen), priorSched(priorSched), alpha(alpha), M_samp(M_samp), B_resamp(B_resamp), w(w), cover(cover) {}   

   // operator
void operator()(std::size_t begin, std::size_t end) {
		for (std::size_t i = begin; i < end; i++) {
			cover[i] = GibbsMCMCkde(nn, data, nnp, priordata, priorweight, thetaboot, bootmean0, bootmean1, bootmeanYI, kdecdflen, kdecdfboot1, kdecdfboot2, kdecdfboot3, kdecdfboot1p, kdecdfboot2p, kdecdfboot3p, scheduleLen, priorSched, alpha, M_samp, w, i);	
		}
	}
};

// [[Rcpp::export]]
NumericVector rcpp_parallel_yi_kde(NumericVector nn, NumericMatrix data, NumericVector nnp, NumericMatrix priordata, NumericVector priorweight, NumericMatrix thetaboot, NumericVector bootmean0,
	NumericVector bootmean1, NumericVector bootmeanYI, NumericVector kdecdflen, NumericMatrix kdecdfboot1, NumericMatrix kdecdfboot2, NumericMatrix kdecdfboot3, NumericMatrix kdecdfboot1p, NumericMatrix kdecdfboot2p, NumericMatrix kdecdfboot3p, NumericVector scheduleLen, NumericMatrix priorSched, NumericVector alpha, NumericVector M_samp, NumericVector B_resamp,
	NumericVector w) {
	
   int B = int(B_resamp[0]);
   // allocate the matrix we will return
   NumericVector cover(B,2.0); 

   // create the worker
   GPCYI_yi_kde_mcmc_parallel gpcWorker(nn, data, nnp, priordata, priorweight, thetaboot, bootmean0, bootmean1, bootmeanYI, kdecdflen, kdecdfboot1, kdecdfboot2, kdecdfboot3, kdecdfboot1p, kdecdfboot2p, kdecdfboot3p, scheduleLen, priorSched, alpha, M_samp, B_resamp, w, cover);
     
   // call it with parallelFor
   
   parallelFor(0, B, gpcWorker);

   return cover;
}


struct GPCYI_yi_mcmc_smooth_parallel : public Worker {
	
	const RVector<double> nn;
	const RMatrix<double> data1;
	const RMatrix<double> data2;
	const RMatrix<double> thetaboot;
	const RVector<double> bootmean0;
	const RVector<double> bootmean1;
	const RMatrix<double> databoot1;
	const RMatrix<double> databoot2;
	const RVector<double> normprior;
	const RVector<double> scheduleLen;
	const RMatrix<double> propSched;
	const RVector<double> ddelta;
	const RVector<double> alpha;
	const RVector<double> M_samp;
	const RVector<double> B_resamp;
	const RVector<double> w;
	RMatrix<double> cover;
	std::vector<double> temp;

   // initialize with source and destination
   GPCYI_yi_mcmc_smooth_parallel(const NumericVector nn, const NumericMatrix data1, const NumericMatrix data2, const NumericMatrix thetaboot,
	const NumericVector bootmean0, const NumericVector bootmean1, const NumericMatrix databoot1, const NumericMatrix databoot2, const NumericVector normprior, const NumericVector scheduleLen, const NumericMatrix propSched,
	const NumericVector ddelta, const NumericVector alpha, const NumericVector M_samp, const NumericVector B_resamp, const NumericVector w, NumericMatrix cover, std::vector<double> temp) 
			: nn(nn), data1(data1), data2(data2), thetaboot(thetaboot), bootmean0(bootmean0), bootmean1(bootmean1), databoot1(databoot1), databoot2(databoot2), normprior(normprior), scheduleLen(scheduleLen), propSched(propSched), ddelta(ddelta), alpha(alpha), M_samp(M_samp), B_resamp(B_resamp), w(w), cover(cover), temp(temp) {}   

   // operator
void operator()(std::size_t begin, std::size_t end) {
		for (std::size_t i = begin; i < end; i++) {
			temp = GibbsMCMCsmooth(nn, data1, data2, thetaboot, bootmean0, bootmean1, databoot1, databoot2, normprior, scheduleLen, propSched, ddelta, alpha, M_samp, w, i);	
			cover(i,0) = temp[0];cover(i,1) = temp[1];
		}
	}
};


// [[Rcpp::export]]
NumericVector rcpp_parallel_smooth_yi(NumericVector nn, NumericMatrix data1, NumericMatrix data2, NumericMatrix thetaboot, NumericVector bootmean0,
	NumericVector bootmean1, NumericMatrix databoot1, NumericMatrix databoot2, NumericVector normprior, NumericVector scheduleLen, NumericMatrix propSched, NumericVector ddelta, NumericVector alpha, NumericVector M_samp, NumericVector B_resamp, 
	NumericVector w) {
		
   int B = int(B_resamp[0]);
   // allocate the matrix we will return
   NumericMatrix cover(B,2); 
   std::vector<double> temp;
   // create the worker
   GPCYI_yi_mcmc_smooth_parallel gpcWorker(nn, data1, data2, thetaboot, bootmean0,
	bootmean1, databoot1, databoot2, normprior, scheduleLen, propSched, ddelta, alpha, M_samp, B_resamp, 
	w, cover, temp);
     
   // call it with parallelFor
   
   parallelFor(0, B, gpcWorker);

   return cover;
}
	
	
				      
// with prior data
				   
struct GPCYI_yi_mcmcp_smooth_parallel : public Worker {
	
	const RVector<double> nn;
	const RMatrix<double> data1;
	const RMatrix<double> data2;
	const RMatrix<double> thetaboot;
	const RVector<double> bootmean0;
	const RVector<double> bootmean1;
	const RMatrix<double> databoot1;
	const RMatrix<double> databoot2;
	const RMatrix<double> priordata1;
	const RMatrix<double> priordata2;
	const RVector<double> priorweight;
	const RVector<double> scheduleLen;
	const RMatrix<double> propSched;
	const RVector<double> ddelta;
	const RVector<double> alpha;
	const RVector<double> M_samp;
	const RVector<double> B_resamp;
	const RVector<double> w;
	RMatrix<double> cover;
	std::vector<double> temp;

   // initialize with source and destination
   GPCYI_yi_mcmcp_smooth_parallel(const NumericVector nn, const NumericMatrix data1, const NumericMatrix data2, const NumericMatrix thetaboot,
	const NumericVector bootmean0, const NumericVector bootmean1, const NumericMatrix databoot1, const NumericMatrix databoot2, const NumericMatrix priordata1, const NumericMatrix priordata2, const NumericVector priorweight, const NumericVector scheduleLen, const NumericMatrix propSched,
	const NumericVector ddelta, const NumericVector alpha, const NumericVector M_samp, const NumericVector B_resamp, const NumericVector w, NumericMatrix cover, std::vector<double> temp) 
			: nn(nn), data1(data1), data2(data2), thetaboot(thetaboot), bootmean0(bootmean0), bootmean1(bootmean1), databoot1(databoot1), databoot2(databoot2), priordata1(priordata1), priordata2(priordata2), priorweight(priorweight), scheduleLen(scheduleLen), propSched(propSched), ddelta(ddelta), alpha(alpha), M_samp(M_samp), B_resamp(B_resamp), w(w), cover(cover), temp(temp) {}   

   // operator
void operator()(std::size_t begin, std::size_t end) {
		for (std::size_t i = begin; i < end; i++) {
			temp = GibbsMCMCpsmooth(nn, data1, data2, thetaboot, bootmean0, bootmean1, databoot1, databoot2, priordata1, priordata2, priorweight, scheduleLen, propSched, ddelta, alpha, M_samp, w, i);	
			cover(i,0) = temp[0];cover(i,1) = temp[1];
		}
	}
};


// [[Rcpp::export]]
NumericVector rcpp_parallel_smoothp_yi(NumericVector nn, NumericMatrix data1, NumericMatrix data2, NumericMatrix thetaboot, NumericVector bootmean0,
	NumericVector bootmean1, NumericMatrix databoot1, NumericMatrix databoot2, NumericMatrix priordata1, NumericMatrix priordata2, NumericVector priorweight, NumericVector scheduleLen, NumericMatrix propSched, NumericVector ddelta, NumericVector alpha, NumericVector M_samp, NumericVector B_resamp, 
	NumericVector w) {
		
   int B = int(B_resamp[0]);
   // allocate the matrix we will return
   NumericMatrix cover(B,2); 
   std::vector<double> temp;
   // create the worker
   GPCYI_yi_mcmcp_smooth_parallel gpcWorker(nn, data1, data2, thetaboot, bootmean0,
	bootmean1, databoot1, databoot2, priordata1, priordata2, priorweight, scheduleLen, propSched, ddelta, alpha, M_samp, B_resamp, 
	w, cover, temp);
     
   // call it with parallelFor
   
   parallelFor(0, B, gpcWorker);

   return cover;
}
	
				      
// with prior data
				   
struct GPCYI_yi_mcmc3_smooth_parallel : public Worker {
	
	const RVector<double> nn;
	const RMatrix<double> data1;
	const RMatrix<double> data2;
	const RMatrix<double> thetaboot;
	const RVector<double> bootmean0;
	const RVector<double> bootmean1;
	const RMatrix<double> databoot1;
	const RMatrix<double> databoot2;
	const RVector<double> priorsched;
	const RVector<double> scheduleLen;
	const RMatrix<double> propSched;
	const RVector<double> ddelta;
	const RVector<double> alpha;
	const RVector<double> M_samp;
	const RVector<double> B_resamp;
	const RVector<double> w;
	RMatrix<double> cover;
	std::vector<double> temp;

   // initialize with source and destination
   GPCYI_yi_mcmc3_smooth_parallel(const NumericVector nn, const NumericMatrix data1, const NumericMatrix data2, const NumericMatrix thetaboot,
	const NumericVector bootmean0, const NumericVector bootmean1, const NumericMatrix databoot1, const NumericMatrix databoot2, const NumericVector priorsched, const NumericVector scheduleLen, const NumericMatrix propSched,
	const NumericVector ddelta, const NumericVector alpha, const NumericVector M_samp, const NumericVector B_resamp, const NumericVector w, NumericMatrix cover, std::vector<double> temp) 
			: nn(nn), data1(data1), data2(data2), thetaboot(thetaboot), bootmean0(bootmean0), bootmean1(bootmean1), databoot1(databoot1), databoot2(databoot2), priorsched(priorsched), scheduleLen(scheduleLen), propSched(propSched), ddelta(ddelta), alpha(alpha), M_samp(M_samp), B_resamp(B_resamp), w(w), cover(cover), temp(temp) {}   

   // operator
void operator()(std::size_t begin, std::size_t end) {
		for (std::size_t i = begin; i < end; i++) {
			temp = GibbsMCMC3loopsmooth(nn, data1, data2, thetaboot, bootmean0, bootmean1, databoot1, databoot2, priorsched, scheduleLen, propSched, ddelta, alpha, M_samp, w, i);	
			cover(i,0) = temp[0];cover(i,1) = temp[1];
		}
	}
};


// [[Rcpp::export]]
NumericVector rcpp_parallel_smooth3_yi(NumericVector nn, NumericMatrix data1, NumericMatrix data2, NumericMatrix thetaboot, NumericVector bootmean0,
	NumericVector bootmean1, NumericMatrix databoot1, NumericMatrix databoot2, NumericVector priorsched, NumericVector scheduleLen, NumericMatrix propSched, NumericVector ddelta, NumericVector alpha, NumericVector M_samp, NumericVector B_resamp, 
	NumericVector w) {
		
   int B = int(B_resamp[0]);
   // allocate the matrix we will return
   NumericMatrix cover(B,2); 
   std::vector<double> temp;
   // create the worker
   GPCYI_yi_mcmc3_smooth_parallel gpcWorker(nn, data1, data2, thetaboot, bootmean0,
	bootmean1, databoot1, databoot2, priorsched, scheduleLen, propSched, ddelta, alpha, M_samp, B_resamp, 
	w, cover, temp);
     
   // call it with parallelFor
   
   parallelFor(0, B, gpcWorker);

   return cover;
}				      
				      
				      
				      
				      

// [[Rcpp::export]]
Rcpp::List GPCYI_yi_parallel(SEXP & nn, SEXP & data, SEXP & nnp, SEXP & priordata, SEXP & priorweight, SEXP & theta_boot, SEXP & data_boot, SEXP & scheduleLen, SEXP & priorSched, SEXP & alpha, SEXP & M_samp, SEXP & B_resamp) {

RNGScope scp;
Rcpp::Function _GPCYI_rcpp_parallel_yi("rcpp_parallel_yi");
List result;
List finalsample;
double eps 			= 0.01; 
NumericVector nn_ = Rcpp::as<NumericVector>(nn);
NumericMatrix data_ = Rcpp::as<NumericMatrix>(data);
NumericVector nnp_ = Rcpp::as<NumericVector>(nnp);
NumericMatrix priordata_ = Rcpp::as<NumericMatrix>(priordata);
NumericVector priorweight_ = Rcpp::as<NumericVector>(priorweight);
NumericMatrix thetaboot_ = Rcpp::as<NumericMatrix>(theta_boot);
NumericVector bootmean0(1,0.0);
NumericVector bootmean1(1,0.0);
NumericMatrix databoot_ = Rcpp::as<NumericMatrix>(data_boot);
NumericVector scheduleLen_ = Rcpp::as<NumericVector>(scheduleLen);
NumericMatrix priorSched_ = Rcpp::as<NumericMatrix>(priorSched);
NumericVector alpha_ = Rcpp::as<NumericVector>(alpha);
NumericVector M_samp_ = Rcpp::as<NumericVector>(M_samp);
NumericVector B_resamp_ = Rcpp::as<NumericVector>(B_resamp);
NumericVector w(1,0.5);
double diff;
bool go 			= TRUE;
int t				=1; 
double sumcover;
int B = int(B_resamp_[0]);
NumericVector cover;
	
for (int i=0; i<B; i++) {
	bootmean0[0] = bootmean0[0] + thetaboot_(i,0);
	bootmean1[0] = bootmean1[0] + thetaboot_(i,1);
}
bootmean0 = bootmean0/B;
bootmean1 = bootmean1/B;

while(go){	
cover = _GPCYI_rcpp_parallel_yi(nn_, data_, nnp_, priordata_, priorweight_, thetaboot_, bootmean0, bootmean1, databoot_, alpha_, M_samp_, B_resamp_, w);
sumcover = 0.0;
for(int s = 0; s<B; s++){sumcover = sumcover + cover(s);}
diff = (sumcover/B) - (1.0-alpha_[0]);
if(((abs(diff)<= eps)&&(diff>=0)) || t>16) {
   go = FALSE;
} else {
   t = t+1;
   w[0] = fmax(w[0] + (pow(1+t,-0.51)*diff),0.1);
} 
}

// Final sample

NumericVector M_final; M_final[0] = 2*M_samp_[0];
finalsample = GibbsMCMC2(nn_, data_, nnp_, priordata_, priorweight_, thetaboot_, bootmean0, bootmean1, scheduleLen_, priorSched_, alpha_, M_final, w);
	
result = Rcpp::List::create(Rcpp::Named("w") = w,Rcpp::Named("t") = t,Rcpp::Named("diff") = diff, Rcpp::Named("list_cis") = finalsample);
	
return result;
}




// [[Rcpp::export]]

inline bool compare(std::array<double, 6> a, std::array<double, 6> b){
    return (a[5] < b[5]);
}





