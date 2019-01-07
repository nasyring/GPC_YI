// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;
using namespace std;



// GPC_yi_parallel
Rcpp::List GPCYI_yi_parallel(SEXP & nn, SEXP & data, SEXP & nnp, SEXP & priordata, SEXP & priorweight, SEXP & theta_boot, SEXP & data_boot, SEXP & alpha, SEXP & M_samp, SEXP & B_resamp);
RcppExport SEXP GPCYI_GPCYI_yi_parallel(SEXP nnSEXP, SEXP dataSEXP, SEXP nnpSEXP, SEXP priordataSEXP, SEXP priorweightSEXP, SEXP theta_bootSEXP, SEXP data_bootSEXP, SEXP alphaSEXP, SEXP M_sampSEXP, SEXP B_resampSEXP){
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
    Rcpp::traits::input_parameter< SEXP & >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< SEXP & >::type M_samp(M_sampSEXP);
    Rcpp::traits::input_parameter< SEXP & >::type B_resamp(B_resampSEXP);
    __result = Rcpp::wrap(GPCYI_yi_parallel(nn, data, nnp, priordata, priorweight, theta_boot, data_boot, alpha, M_samp, B_resamp));
    return __result;
END_RCPP
}





// rcpp_parallel_yi
NumericVector rcpp_parallel_yi(NumericVector nn, NumericMatrix data, NumericVector nnp, NumericMatrix priordata, NumericVector priorweight, NumericMatrix thetaboot, NumericVector bootmean0,
	NumericVector bootmean1, NumericMatrix databoot, NumericVector alpha, NumericVector M_samp, NumericVector B_resamp,
	NumericVector w);
RcppExport SEXP GPCYI_rcpp_parallel_yi(SEXP nnSEXP, SEXP dataSEXP, SEXP nnpSEXP, SEXP priordataSEXP, SEXP priorweightSEXP, SEXP thetabootSEXP, SEXP bootmean0SEXP, SEXP bootmean1SEXP, SEXP databootSEXP,
                                    SEXP alphaSEXP, SEXP M_sampSEXP, SEXP B_resampSEXP, SEXP wSEXP){
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
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type M_samp(M_sampSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type B_resamp(B_resampSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    __result = Rcpp::wrap(rcpp_parallel_yi(nn, data, nnp, priordata, priorweight, thetaboot, bootmean0, bootmean1, databoot, alpha, M_samp, B_resamp, w));
    return __result;
END_RCPP
}    



Rcpp::List GibbsMCMC2(NumericVector nn, NumericMatrix data, NumericVector nnp, NumericMatrix priordata, NumericVector priorweight, NumericMatrix thetaboot,
	NumericVector bootmean0, NumericVector bootmean1, NumericVector alpha, NumericVector M_samp, NumericVector w);
RcppExport SEXP GPCYI_GibbsMCMC2(SEXP nnSEXP, SEXP dataSEXP, SEXP nnpSEXP, SEXP priordataSEXP, SEXP priorweightSEXP, SEXP thetabootSEXP, SEXP bootmean0SEXP, SEXP bootmean1SEXP, SEXP alphaSEXP, SEXP M_sampSEXP, SEXP wSEXP){
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
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type M_samp(M_sampSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    __result = Rcpp::wrap(GibbsMCMC2(nn, data, nnp, priordata, priorweight, thetaboot, bootmean0, bootmean1, alpha, M_samp, w));
    return __result;
END_RCPP
}




static const R_CallMethodDef CallEntries[] = {
    {"GPCYI_GPCYI_yi_parallel", (DL_FUNC) &GPCYI_GPCYI_yi_parallel, 10},
    {"GPCYI_rcpp_parallel_yi", (DL_FUNC) &GPCYI_rcpp_parallel_yi, 13},
    {"GPCYI_GibbsMCMC2", (DL_FUNC) &GPCYI_GibbsMCMC2, 11},
    {NULL, NULL, 0}
};

RcppExport void R_init_GPCYI(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
    



