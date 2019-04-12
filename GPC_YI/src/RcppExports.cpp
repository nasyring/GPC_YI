// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;
using namespace std;



// GPC_yi_parallel
Rcpp::List GPCYI_yi_parallel(SEXP & nn, SEXP & data, SEXP & nnp, SEXP & priordata, SEXP & priorweight, SEXP & theta_boot, SEXP & data_boot, SEXP & scheduleLen, SEXP & priorSched, SEXP & alpha, SEXP & M_samp, SEXP & B_resamp);
RcppExport SEXP GPCYI_GPCYI_yi_parallel(SEXP nnSEXP, SEXP dataSEXP, SEXP nnpSEXP, SEXP priordataSEXP, SEXP priorweightSEXP, SEXP theta_bootSEXP, SEXP data_bootSEXP, SEXP scheduleLenSEXP, SEXP priorSchedSEXP, SEXP alphaSEXP, SEXP M_sampSEXP, SEXP B_resampSEXP){
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< SEXP & >::type nn(nnSEXP);
    Rcpp::traits::input_parameter< SEXP & >::type data(dataSEXP);
    Rcpp::traits::input_parameter< SEXP & >::type nnp(nnpSEXP);
    Rcpp::traits::input_parameter< SEXP & >::type priordata(priordataSEXP);
    Rcpp::traits::input_parameter< SEXP & >::type priorweight(priorweightSEXP);
    Rcpp::traits::input_parameter< SEXP & >::type theta_boot(theta_bootSEXP);
    Rcpp::traits::input_parameter< SEXP & >::type data_boot(data_bootSEXP);
    Rcpp::traits::input_parameter< SEXP & >::type scheduleLen(scheduleLenSEXP);
    Rcpp::traits::input_parameter< SEXP & >::type priorSched(priorSchedSEXP);
    Rcpp::traits::input_parameter< SEXP & >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< SEXP & >::type M_samp(M_sampSEXP);
    Rcpp::traits::input_parameter< SEXP & >::type B_resamp(B_resampSEXP);
    __result = Rcpp::wrap(GPCYI_yi_parallel(nn, data, nnp, priordata, priorweight, theta_boot, data_boot, scheduleLen, priorSched, alpha, M_samp, B_resamp));
    return __result;
END_RCPP
}


   	
Rcpp::List GridSearchKDE(int n, NumericVector mesh1, NumericVector mesh2, NumericVector mesh3, NumericVector cdf1, NumericVector cdf2, NumericVector cdf3, int n12, int n23, NumericVector d12, NumericVector d23 );
RcppExport SEXP GPCYI_GridSearchKDE(SEXP nSEXP, SEXP mesh1SEXP, SEXP mesh2SEXP, SEXP mesh3SEXP, SEXP cdf1SEXP, SEXP cdf2SEXP, SEXP cdf3SEXP, SEXP n12SEXP, SEXP n23SEXP, SEXP d12SEXP, SEXP d23SEXP ){
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mesh1(mesh1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mesh2(mesh2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mesh3(mesh3SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cdf1(cdf1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cdf2(cdf2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cdf3(cdf3SEXP);
    Rcpp::traits::input_parameter< int >::type n12(n12SEXP);
    Rcpp::traits::input_parameter< int >::type n23(n23SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type d12(d12SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type d23(d23SEXP);
    __result = Rcpp::wrap(GridSearchKDE(n, mesh1, mesh2, mesh3, cdf1, cdf2, cdf3, n12, n23, d12, d23));
    return __result;
END_RCPP
}

// rcpp_parallel_yi
NumericVector rcpp_parallel_yi(NumericVector nn, NumericMatrix data, NumericVector nnp, NumericMatrix priordata, NumericVector priorweight, NumericMatrix thetaboot, NumericVector bootmean0,
	NumericVector bootmean1, NumericMatrix databoot, NumericVector scheduleLen, NumericMatrix priorSched, NumericVector alpha, NumericVector M_samp, NumericVector B_resamp,
	NumericVector w);
RcppExport SEXP GPCYI_rcpp_parallel_yi(SEXP nnSEXP, SEXP dataSEXP, SEXP nnpSEXP, SEXP priordataSEXP, SEXP priorweightSEXP, SEXP thetabootSEXP, SEXP bootmean0SEXP, SEXP bootmean1SEXP, SEXP databootSEXP,
                                    SEXP scheduleLenSEXP, SEXP priorSchedSEXP, SEXP alphaSEXP, SEXP M_sampSEXP, SEXP B_resampSEXP, SEXP wSEXP){
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type nn(nnSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nnp(nnpSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type priordata(priordataSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type priorweight(priorweightSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type thetaboot(thetabootSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bootmean0(bootmean0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bootmean1(bootmean1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type databoot(databootSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type scheduleLen(scheduleLenSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type priorSched(priorSchedSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type M_samp(M_sampSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type B_resamp(B_resampSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    __result = Rcpp::wrap(rcpp_parallel_yi(nn, data, nnp, priordata, priorweight, thetaboot, bootmean0, bootmean1, databoot, scheduleLen, priorSched, alpha, M_samp, B_resamp, w));
    return __result;
END_RCPP
}    


// rcpp_parallel_yi
NumericVector rcpp_parallel_smooth_yi(NumericVector nn, NumericMatrix data1, NumericMatrix data2, NumericMatrix thetaboot, NumericVector bootmean0,
	NumericVector bootmean1, NumericMatrix databoot1, NumericMatrix databoot2, NumericVector normprior, NumericVector scheduleLen, NumericMatrix propSched, NumericVector ddelta, NumericVector alpha, NumericVector M_samp, NumericVector B_resamp, 
	NumericVector w);
RcppExport SEXP GPCYI_rcpp_parallel_smooth_yi(SEXP nnSEXP, SEXP data1SEXP, SEXP data2SEXP, SEXP thetabootSEXP, SEXP bootmean0SEXP,
	SEXP bootmean1SEXP, SEXP databoot1SEXP, SEXP databoot2SEXP, SEXP normpriorSEXP, SEXP scheduleLenSEXP, SEXP propSchedSEXP, SEXP ddeltaSEXP, SEXP alphaSEXP, SEXP M_sampSEXP, SEXP B_resampSEXP, 
	SEXP wSEXP){
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type nn(nnSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type data1(data1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type data2(data2SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type thetaboot(thetabootSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bootmean0(bootmean0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bootmean1(bootmean1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type databoot1(databoot1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type databoot2(databoot2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type normprior(normpriorSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type scheduleLen(scheduleLenSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type propSched(propSchedSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ddelta(ddeltaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type M_samp(M_sampSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type B_resamp(B_resampSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    __result = Rcpp::wrap(rcpp_parallel_smooth_yi(nn, data1, data2, thetaboot, bootmean0, bootmean1, databoot1, databoot2, normprior, scheduleLen, propSched, ddelta, alpha, M_samp, B_resamp, w));
    return __result;
END_RCPP
}    

// rcpp_parallel_yi smooth prior data
NumericVector rcpp_parallel_smoothp_yi(NumericVector nn, NumericMatrix data1, NumericMatrix data2, NumericMatrix thetaboot, NumericVector bootmean0,
	NumericVector bootmean1, NumericMatrix databoot1, NumericMatrix databoot2, NumericMatrix priordata1, NumericMatrix priordata2, NumericVector priorweight, NumericVector scheduleLen, NumericMatrix propSched, NumericVector ddelta, NumericVector alpha, NumericVector M_samp, NumericVector B_resamp, 
	NumericVector w);
RcppExport SEXP GPCYI_rcpp_parallel_smoothp_yi(SEXP nnSEXP, SEXP data1SEXP, SEXP data2SEXP, SEXP thetabootSEXP, SEXP bootmean0SEXP,
	SEXP bootmean1SEXP, SEXP databoot1SEXP, SEXP databoot2SEXP, SEXP priordata1SEXP, SEXP priordata2SEXP, SEXP priorweightSEXP, SEXP scheduleLenSEXP, SEXP propSchedSEXP, SEXP ddeltaSEXP, SEXP alphaSEXP, SEXP M_sampSEXP, SEXP B_resampSEXP, 
	SEXP wSEXP){
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type nn(nnSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type data1(data1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type data2(data2SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type thetaboot(thetabootSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bootmean0(bootmean0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bootmean1(bootmean1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type databoot1(databoot1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type databoot2(databoot2SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type priordata1(priordata1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type priordata2(priordata2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type priorweight(priorweightSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type scheduleLen(scheduleLenSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type propSched(propSchedSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ddelta(ddeltaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type M_samp(M_sampSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type B_resamp(B_resampSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    __result = Rcpp::wrap(rcpp_parallel_smoothp_yi(nn, data1, data2, thetaboot, bootmean0, bootmean1, databoot1, databoot2, priordata1, priordata2, priorweight, scheduleLen, propSched, ddelta, alpha, M_samp, B_resamp, w));
    return __result;
END_RCPP
}  


// rcpp_parallel_yi smooth prior data
NumericVector rcpp_parallel_smooth3_yi(NumericVector nn, NumericMatrix data1, NumericMatrix data2, NumericMatrix thetaboot, NumericVector bootmean0,
	NumericVector bootmean1, NumericMatrix databoot1, NumericMatrix databoot2, NumericVector priorsched, NumericVector scheduleLen, NumericMatrix propSched, NumericVector ddelta, NumericVector alpha, NumericVector M_samp, NumericVector B_resamp, 
	NumericVector w);
RcppExport SEXP GPCYI_rcpp_parallel_smooth3_yi(SEXP nnSEXP, SEXP data1SEXP, SEXP data2SEXP, SEXP thetabootSEXP, SEXP bootmean0SEXP,
	SEXP bootmean1SEXP, SEXP databoot1SEXP, SEXP databoot2SEXP, SEXP priorschedSEXP, SEXP scheduleLenSEXP, SEXP propSchedSEXP, SEXP ddeltaSEXP, SEXP alphaSEXP, SEXP M_sampSEXP, SEXP B_resampSEXP, 
	SEXP wSEXP){
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type nn(nnSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type data1(data1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type data2(data2SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type thetaboot(thetabootSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bootmean0(bootmean0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bootmean1(bootmean1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type databoot1(databoot1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type databoot2(databoot2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type priorsched(priorschedSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type scheduleLen(scheduleLenSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type propSched(propSchedSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ddelta(ddeltaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type M_samp(M_sampSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type B_resamp(B_resampSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    __result = Rcpp::wrap(rcpp_parallel_smooth3_yi(nn, data1, data2, thetaboot, bootmean0, bootmean1, databoot1, databoot2, priorsched, scheduleLen, propSched, ddelta, alpha, M_samp, B_resamp, w));
    return __result;
END_RCPP
}  

// rcpp_parallel_yi
NumericVector rcpp_parallel_yi_kde(NumericVector nn, NumericMatrix data, NumericVector nnp, NumericMatrix priordata, NumericVector priorweight, NumericMatrix thetaboot, NumericVector bootmean0,
	NumericVector bootmean1, NumericVector bootmeanYI, NumericVector kdecdflen, NumericMatrix kdecdfboot1, NumericMatrix kdecdfboot2, NumericMatrix kdecdfboot3, NumericMatrix kdecdfboot1p, NumericMatrix kdecdfboot2p, NumericMatrix kdecdfboot3p, NumericVector scheduleLen, NumericMatrix priorSched, NumericVector alpha, NumericVector M_samp, NumericVector B_resamp,
	NumericVector w);
RcppExport SEXP GPCYI_rcpp_parallel_yi_kde(SEXP nnSEXP, SEXP dataSEXP, SEXP nnpSEXP, SEXP priordataSEXP, SEXP priorweightSEXP, SEXP thetabootSEXP, SEXP bootmean0SEXP, SEXP bootmean1SEXP, SEXP bootmeanYISEXP, SEXP kdecdflenSEXP,
					   SEXP kdecdfboot1SEXP, SEXP kdecdfboot2SEXP, SEXP kdecdfboot3SEXP, SEXP kdecdfboot1pSEXP, SEXP kdecdfboot2pSEXP, SEXP kdecdfboot3pSEXP, 
                                    SEXP scheduleLenSEXP, SEXP priorSchedSEXP, SEXP alphaSEXP, SEXP M_sampSEXP, SEXP B_resampSEXP, SEXP wSEXP){
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type nn(nnSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nnp(nnpSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type priordata(priordataSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type priorweight(priorweightSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type thetaboot(thetabootSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bootmean0(bootmean0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bootmean1(bootmean1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bootmeanYI(bootmeanYISEXP);
    Rcpp::traits::input_parameter< NumericVector >::type kdecdflen(kdecdflenSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type kdecdfboot1(kdecdfboot1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type kdecdfboot2(kdecdfboot2SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type kdecdfboot3(kdecdfboot3SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type kdecdfboot1p(kdecdfboot1pSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type kdecdfboot2p(kdecdfboot2pSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type kdecdfboot3p(kdecdfboot3pSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type scheduleLen(scheduleLenSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type priorSched(priorSchedSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type M_samp(M_sampSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type B_resamp(B_resampSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    __result = Rcpp::wrap(rcpp_parallel_yi_kde(nn, data, nnp, priordata, priorweight, thetaboot, bootmean0, bootmean1, bootmeanYI, kdecdflen, kdecdfboot1, kdecdfboot2, kdecdfboot3, kdecdfboot1p, kdecdfboot2p, kdecdfboot3p, scheduleLen, priorSched, alpha, M_samp, B_resamp, w));
    return __result;
END_RCPP
}    

			
		
Rcpp::List GibbsMCMC2smooth(NumericVector nn, NumericMatrix data1, NumericMatrix data2, NumericVector bootmean0, NumericVector bootmean1, NumericVector normprior, NumericVector scheduleLen, NumericMatrix propSched, NumericVector alpha, NumericVector ddelta, NumericVector M_samp, NumericVector w);
RcppExport SEXP GPCYI_GibbsMCMC2smooth(SEXP nnSEXP, SEXP data1SEXP, SEXP data2SEXP, SEXP bootmean0SEXP, SEXP bootmean1SEXP, SEXP normpriorSEXP, SEXP scheduleLenSEXP, SEXP propSchedSEXP, SEXP alphaSEXP, SEXP ddeltaSEXP, SEXP M_sampSEXP, SEXP wSEXP){
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type nn(nnSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type data1(data1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type data2(data2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bootmean0(bootmean0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bootmean1(bootmean1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type normprior(normpriorSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type scheduleLen(scheduleLenSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type propSched(propSchedSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ddelta(ddeltaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type M_samp(M_sampSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    __result = Rcpp::wrap(GibbsMCMC2smooth(nn, data1, data2, bootmean0, bootmean1, normprior, scheduleLen, propSched, alpha, ddelta, M_samp, w));
    return __result;
END_RCPP
}

		
Rcpp::List GibbsMCMC3smooth(NumericVector nn, NumericMatrix data1, NumericMatrix data2, NumericVector bootmean0, NumericVector bootmean1, NumericVector priorsched, NumericVector scheduleLen, NumericMatrix propSched, NumericVector alpha, NumericVector ddelta, NumericVector M_samp, NumericVector w);
RcppExport SEXP GPCYI_GibbsMCMC3smooth(SEXP nnSEXP, SEXP data1SEXP, SEXP data2SEXP, SEXP bootmean0SEXP, SEXP bootmean1SEXP, SEXP priorschedSEXP, SEXP scheduleLenSEXP, SEXP propSchedSEXP, SEXP alphaSEXP, SEXP ddeltaSEXP, SEXP M_sampSEXP, SEXP wSEXP){
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type nn(nnSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type data1(data1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type data2(data2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bootmean0(bootmean0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bootmean1(bootmean1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type priorsched(priorschedSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type scheduleLen(scheduleLenSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type propSched(propSchedSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ddelta(ddeltaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type M_samp(M_sampSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    __result = Rcpp::wrap(GibbsMCMC3smooth(nn, data1, data2, bootmean0, bootmean1, priorsched, scheduleLen, propSched, alpha, ddelta, M_samp, w));
    return __result;
END_RCPP
}

Rcpp::List GibbsMCMCp2smooth(NumericVector nn, NumericMatrix data1, NumericMatrix data2, NumericVector bootmean0, NumericVector bootmean1, NumericMatrix priordata1, NumericMatrix priordata2, NumericVector priorweight, NumericVector scheduleLen, NumericMatrix propSched, NumericVector alpha, NumericVector ddelta, NumericVector M_samp, NumericVector w);
RcppExport SEXP GPCYI_GibbsMCMCp2smooth(SEXP nnSEXP, SEXP data1SEXP, SEXP data2SEXP, SEXP bootmean0SEXP, SEXP bootmean1SEXP, SEXP priordata1SEXP, SEXP priordata2SEXP, SEXP priorweightSEXP, SEXP scheduleLenSEXP, SEXP propSchedSEXP, SEXP alphaSEXP, SEXP ddeltaSEXP, SEXP M_sampSEXP, SEXP wSEXP){
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type nn(nnSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type data1(data1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type data2(data2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bootmean0(bootmean0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bootmean1(bootmean1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type priordata1(priordata1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type priordata2(priordata2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type priorweight(priorweightSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type scheduleLen(scheduleLenSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type propSched(propSchedSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ddelta(ddeltaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type M_samp(M_sampSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    __result = Rcpp::wrap(GibbsMCMCp2smooth(nn, data1, data2, bootmean0, bootmean1, priordata1, priordata2, priorweight, scheduleLen, propSched, alpha, ddelta, M_samp, w));
    return __result;
END_RCPP
}
		
		

Rcpp::List GibbsMCMC2(NumericVector nn, NumericMatrix data, NumericVector nnp, NumericMatrix priordata, NumericVector priorweight, NumericMatrix thetaboot,
	NumericVector bootmean0, NumericVector bootmean1, NumericVector scheduleLen, NumericMatrix priorSched, NumericVector alpha, NumericVector M_samp, NumericVector w);
RcppExport SEXP GPCYI_GibbsMCMC2(SEXP nnSEXP, SEXP dataSEXP, SEXP nnpSEXP, SEXP priordataSEXP, SEXP priorweightSEXP, SEXP thetabootSEXP, SEXP bootmean0SEXP, SEXP bootmean1SEXP, SEXP scheduleLenSEXP, SEXP priorSchedSEXP, SEXP alphaSEXP, SEXP M_sampSEXP, SEXP wSEXP){
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type nn(nnSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nnp(nnpSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type priordata(priordataSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type priorweight(priorweightSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type thetaboot(thetabootSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bootmean0(bootmean0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bootmean1(bootmean1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type scheduleLen(scheduleLenSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type priorSched(priorSchedSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type M_samp(M_sampSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    __result = Rcpp::wrap(GibbsMCMC2(nn, data, nnp, priordata, priorweight, thetaboot, bootmean0, bootmean1, scheduleLen, priorSched, alpha, M_samp, w));
    return __result;
END_RCPP
}

Rcpp::List GibbsMCMCkde2(NumericVector nn, NumericMatrix data, NumericVector nnp, NumericMatrix priordata, NumericVector priorweight, NumericMatrix thetaboot,
	NumericVector bootmean0, NumericVector bootmean1, NumericVector kdecdflen,  NumericMatrix kdecdfboot1, NumericMatrix kdecdfboot2, NumericMatrix kdecdfboot3, NumericMatrix kdecdfboot1p, NumericMatrix kdecdfboot2p, NumericMatrix kdecdfboot3p, NumericVector scheduleLen, NumericMatrix priorSched, NumericVector alpha, NumericVector M_samp, NumericVector w);
RcppExport SEXP GPCYI_GibbsMCMCkde2(SEXP nnSEXP, SEXP dataSEXP, SEXP nnpSEXP, SEXP priordataSEXP, SEXP priorweightSEXP, SEXP thetabootSEXP, SEXP bootmean0SEXP, SEXP bootmean1SEXP, SEXP kdecdflenSEXP, SEXP kdecdfboot1SEXP, SEXP kdecdfboot2SEXP, SEXP kdecdfboot3SEXP, SEXP kdecdfboot1pSEXP, SEXP kdecdfboot2pSEXP, SEXP kdecdfboot3pSEXP,
				    SEXP scheduleLenSEXP, SEXP priorSchedSEXP, SEXP alphaSEXP, SEXP M_sampSEXP, SEXP wSEXP){
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type nn(nnSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nnp(nnpSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type priordata(priordataSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type priorweight(priorweightSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type thetaboot(thetabootSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bootmean0(bootmean0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bootmean1(bootmean1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type kdecdflen(kdecdflenSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type kdecdfboot1(kdecdfboot1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type kdecdfboot2(kdecdfboot2SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type kdecdfboot3(kdecdfboot3SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type kdecdfboot1p(kdecdfboot1pSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type kdecdfboot2p(kdecdfboot2pSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type kdecdfboot3p(kdecdfboot3pSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type scheduleLen(scheduleLenSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type priorSched(priorSchedSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type M_samp(M_sampSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    __result = Rcpp::wrap(GibbsMCMCkde2(nn, data, nnp, priordata, priorweight, thetaboot, bootmean0, bootmean1, kdecdflen, kdecdfboot1, kdecdfboot2, kdecdfboot3, kdecdfboot1p, kdecdfboot2p, kdecdfboot3p, scheduleLen, priorSched, alpha, M_samp, w));
    return __result;
END_RCPP
}


static const R_CallMethodDef CallEntries[] = {
    {"GPCYI_GPCYI_yi_parallel", (DL_FUNC) &GPCYI_GPCYI_yi_parallel, 12},
    {"GPCYI_GridSearchKDE", (DL_FUNC) &GPCYI_GridSearchKDE, 11},
    {"GPCYI_rcpp_parallel_yi", (DL_FUNC) &GPCYI_rcpp_parallel_yi, 15},
    {"GPCYI_rcpp_parallel_smooth_yi", (DL_FUNC) &GPCYI_rcpp_parallel_smooth_yi, 16},
    {"GPCYI_rcpp_parallel_smooth3_yi", (DL_FUNC) &GPCYI_rcpp_parallel_smooth3_yi, 16},
    {"GPCYI_rcpp_parallel_smoothp_yi", (DL_FUNC) &GPCYI_rcpp_parallel_smoothp_yi, 18},
    {"GPCYI_rcpp_parallel_yi_kde", (DL_FUNC) &GPCYI_rcpp_parallel_yi_kde, 22},
    {"GPCYI_GibbsMCMC2", (DL_FUNC) &GPCYI_GibbsMCMC2, 13},
    {"GPCYI_GibbsMCMC2smooth", (DL_FUNC) &GPCYI_GibbsMCMC2smooth, 12},
    {"GPCYI_GibbsMCMC3smooth", (DL_FUNC) &GPCYI_GibbsMCMC2smooth, 12},
    {"GPCYI_GibbsMCMCp2smooth", (DL_FUNC) &GPCYI_GibbsMCMCp2smooth, 14},
    {"GPCYI_GibbsMCMCkde2", (DL_FUNC) &GPCYI_GibbsMCMCkde2, 20},
    {NULL, NULL, 0}
};

RcppExport void R_init_GPCYI(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
    



