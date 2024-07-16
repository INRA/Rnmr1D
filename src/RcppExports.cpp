// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// lorentz
double lorentz(double x, double x0, double s, double a);
RcppExport SEXP _Rnmr1D_lorentz(SEXP xSEXP, SEXP x0SEXP, SEXP sSEXP, SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(lorentz(x, x0, s, a));
    return rcpp_result_gen;
END_RCPP
}
// gauss
double gauss(double x, double x0, double s, double a);
RcppExport SEXP _Rnmr1D_gauss(SEXP xSEXP, SEXP x0SEXP, SEXP sSEXP, SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(gauss(x, x0, s, a));
    return rcpp_result_gen;
END_RCPP
}
// pvoigt
double pvoigt(double x, double x0, double s, double a, double eta);
RcppExport SEXP _Rnmr1D_pvoigt(SEXP xSEXP, SEXP x0SEXP, SEXP sSEXP, SEXP aSEXP, SEXP etaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type eta(etaSEXP);
    rcpp_result_gen = Rcpp::wrap(pvoigt(x, x0, s, a, eta));
    return rcpp_result_gen;
END_RCPP
}
// C_fSavGol
SEXP C_fSavGol(SEXP s, int m, int nl, int nr);
RcppExport SEXP _Rnmr1D_C_fSavGol(SEXP sSEXP, SEXP mSEXP, SEXP nlSEXP, SEXP nrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type s(sSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type nl(nlSEXP);
    Rcpp::traits::input_parameter< int >::type nr(nrSEXP);
    rcpp_result_gen = Rcpp::wrap(C_fSavGol(s, m, nl, nr));
    return rcpp_result_gen;
END_RCPP
}
// C_FilterbyWT
SEXP C_FilterbyWT(SEXP s, int wavelet, double threshold, int verbose);
RcppExport SEXP _Rnmr1D_C_FilterbyWT(SEXP sSEXP, SEXP waveletSEXP, SEXP thresholdSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type s(sSEXP);
    Rcpp::traits::input_parameter< int >::type wavelet(waveletSEXP);
    Rcpp::traits::input_parameter< double >::type threshold(thresholdSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(C_FilterbyWT(s, wavelet, threshold, verbose));
    return rcpp_result_gen;
END_RCPP
}
// C_FilterbyThreshold
SEXP C_FilterbyThreshold(SEXP s, int wavelet, int threshold, int verbose);
RcppExport SEXP _Rnmr1D_C_FilterbyThreshold(SEXP sSEXP, SEXP waveletSEXP, SEXP thresholdSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type s(sSEXP);
    Rcpp::traits::input_parameter< int >::type wavelet(waveletSEXP);
    Rcpp::traits::input_parameter< int >::type threshold(thresholdSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(C_FilterbyThreshold(s, wavelet, threshold, verbose));
    return rcpp_result_gen;
END_RCPP
}
// C_Lorentz
SEXP C_Lorentz(SEXP ppm, double amp, double x0, double sigma, double a);
RcppExport SEXP _Rnmr1D_C_Lorentz(SEXP ppmSEXP, SEXP ampSEXP, SEXP x0SEXP, SEXP sigmaSEXP, SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type ppm(ppmSEXP);
    Rcpp::traits::input_parameter< double >::type amp(ampSEXP);
    Rcpp::traits::input_parameter< double >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(C_Lorentz(ppm, amp, x0, sigma, a));
    return rcpp_result_gen;
END_RCPP
}
// C_PVoigt
SEXP C_PVoigt(SEXP ppm, double amp, double x0, double sigma, double a, double eta);
RcppExport SEXP _Rnmr1D_C_PVoigt(SEXP ppmSEXP, SEXP ampSEXP, SEXP x0SEXP, SEXP sigmaSEXP, SEXP aSEXP, SEXP etaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type ppm(ppmSEXP);
    Rcpp::traits::input_parameter< double >::type amp(ampSEXP);
    Rcpp::traits::input_parameter< double >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type eta(etaSEXP);
    rcpp_result_gen = Rcpp::wrap(C_PVoigt(ppm, amp, x0, sigma, a, eta));
    return rcpp_result_gen;
END_RCPP
}
// C_OneVoigt
SEXP C_OneVoigt(SEXP X, SEXP Y, SEXP par);
RcppExport SEXP _Rnmr1D_C_OneVoigt(SEXP XSEXP, SEXP YSEXP, SEXP parSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type X(XSEXP);
    Rcpp::traits::input_parameter< SEXP >::type Y(YSEXP);
    Rcpp::traits::input_parameter< SEXP >::type par(parSEXP);
    rcpp_result_gen = Rcpp::wrap(C_OneVoigt(X, Y, par));
    return rcpp_result_gen;
END_RCPP
}
// C_peakFinder
SEXP C_peakFinder(SEXP spec, SEXP ppmrange, Nullable<List> filt, Nullable<List> params, int verbose);
RcppExport SEXP _Rnmr1D_C_peakFinder(SEXP specSEXP, SEXP ppmrangeSEXP, SEXP filtSEXP, SEXP paramsSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type spec(specSEXP);
    Rcpp::traits::input_parameter< SEXP >::type ppmrange(ppmrangeSEXP);
    Rcpp::traits::input_parameter< Nullable<List> >::type filt(filtSEXP);
    Rcpp::traits::input_parameter< Nullable<List> >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(C_peakFinder(spec, ppmrange, filt, params, verbose));
    return rcpp_result_gen;
END_RCPP
}
// C_peakOptimize
SEXP C_peakOptimize(SEXP spec, SEXP ppmrange, SEXP params, int verbose);
RcppExport SEXP _Rnmr1D_C_peakOptimize(SEXP specSEXP, SEXP ppmrangeSEXP, SEXP paramsSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type spec(specSEXP);
    Rcpp::traits::input_parameter< SEXP >::type ppmrange(ppmrangeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(C_peakOptimize(spec, ppmrange, params, verbose));
    return rcpp_result_gen;
END_RCPP
}
// C_peakFiltering
SEXP C_peakFiltering(SEXP spec, SEXP peaks, double ratioPN);
RcppExport SEXP _Rnmr1D_C_peakFiltering(SEXP specSEXP, SEXP peaksSEXP, SEXP ratioPNSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type spec(specSEXP);
    Rcpp::traits::input_parameter< SEXP >::type peaks(peaksSEXP);
    Rcpp::traits::input_parameter< double >::type ratioPN(ratioPNSEXP);
    rcpp_result_gen = Rcpp::wrap(C_peakFiltering(spec, peaks, ratioPN));
    return rcpp_result_gen;
END_RCPP
}
// C_specModel
SEXP C_specModel(SEXP spec, SEXP ppmrange, SEXP peaks);
RcppExport SEXP _Rnmr1D_C_specModel(SEXP specSEXP, SEXP ppmrangeSEXP, SEXP peaksSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type spec(specSEXP);
    Rcpp::traits::input_parameter< SEXP >::type ppmrange(ppmrangeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type peaks(peaksSEXP);
    rcpp_result_gen = Rcpp::wrap(C_specModel(spec, ppmrange, peaks));
    return rcpp_result_gen;
END_RCPP
}
// C_MyFuncTest2
SEXP C_MyFuncTest2(SEXP spec, int n1, int n2);
RcppExport SEXP _Rnmr1D_C_MyFuncTest2(SEXP specSEXP, SEXP n1SEXP, SEXP n2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type spec(specSEXP);
    Rcpp::traits::input_parameter< int >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< int >::type n2(n2SEXP);
    rcpp_result_gen = Rcpp::wrap(C_MyFuncTest2(spec, n1, n2));
    return rcpp_result_gen;
END_RCPP
}
// SDL
SEXP SDL(SEXP x, double Sigma);
RcppExport SEXP _Rnmr1D_SDL(SEXP xSEXP, SEXP SigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type Sigma(SigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(SDL(x, Sigma));
    return rcpp_result_gen;
END_RCPP
}
// C_write_pack
void C_write_pack(SEXP x, double pmin, double pmax, SEXP ff);
RcppExport SEXP _Rnmr1D_C_write_pack(SEXP xSEXP, SEXP pminSEXP, SEXP pmaxSEXP, SEXP ffSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type pmin(pminSEXP);
    Rcpp::traits::input_parameter< double >::type pmax(pmaxSEXP);
    Rcpp::traits::input_parameter< SEXP >::type ff(ffSEXP);
    C_write_pack(x, pmin, pmax, ff);
    return R_NilValue;
END_RCPP
}
// C_read_pack
SEXP C_read_pack(SEXP ff);
RcppExport SEXP _Rnmr1D_C_read_pack(SEXP ffSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type ff(ffSEXP);
    rcpp_result_gen = Rcpp::wrap(C_read_pack(ff));
    return rcpp_result_gen;
END_RCPP
}
// C_GlobSeg
SEXP C_GlobSeg(SEXP v, int dN, double sig);
RcppExport SEXP _Rnmr1D_C_GlobSeg(SEXP vSEXP, SEXP dNSEXP, SEXP sigSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type dN(dNSEXP);
    Rcpp::traits::input_parameter< double >::type sig(sigSEXP);
    rcpp_result_gen = Rcpp::wrap(C_GlobSeg(v, dN, sig));
    return rcpp_result_gen;
END_RCPP
}
// lowpass1
SEXP lowpass1(SEXP x, double alpha);
RcppExport SEXP _Rnmr1D_lowpass1(SEXP xSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(lowpass1(x, alpha));
    return rcpp_result_gen;
END_RCPP
}
// WinMoy
double WinMoy(SEXP v, int n1, int n2);
RcppExport SEXP _Rnmr1D_WinMoy(SEXP vSEXP, SEXP n1SEXP, SEXP n2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< int >::type n2(n2SEXP);
    rcpp_result_gen = Rcpp::wrap(WinMoy(v, n1, n2));
    return rcpp_result_gen;
END_RCPP
}
// Smooth
SEXP Smooth(SEXP v, int n);
RcppExport SEXP _Rnmr1D_Smooth(SEXP vSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(Smooth(v, n));
    return rcpp_result_gen;
END_RCPP
}
// fitLines
void fitLines(SEXP s, SEXP b, int n1, int n2);
RcppExport SEXP _Rnmr1D_fitLines(SEXP sSEXP, SEXP bSEXP, SEXP n1SEXP, SEXP n2SEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type s(sSEXP);
    Rcpp::traits::input_parameter< SEXP >::type b(bSEXP);
    Rcpp::traits::input_parameter< int >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< int >::type n2(n2SEXP);
    fitLines(s, b, n1, n2);
    return R_NilValue;
END_RCPP
}
// C_Estime_LB
SEXP C_Estime_LB(SEXP s, int istart, int iend, double WS, double NEIGH, double sig);
RcppExport SEXP _Rnmr1D_C_Estime_LB(SEXP sSEXP, SEXP istartSEXP, SEXP iendSEXP, SEXP WSSEXP, SEXP NEIGHSEXP, SEXP sigSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type s(sSEXP);
    Rcpp::traits::input_parameter< int >::type istart(istartSEXP);
    Rcpp::traits::input_parameter< int >::type iend(iendSEXP);
    Rcpp::traits::input_parameter< double >::type WS(WSSEXP);
    Rcpp::traits::input_parameter< double >::type NEIGH(NEIGHSEXP);
    Rcpp::traits::input_parameter< double >::type sig(sigSEXP);
    rcpp_result_gen = Rcpp::wrap(C_Estime_LB(s, istart, iend, WS, NEIGH, sig));
    return rcpp_result_gen;
END_RCPP
}
// C_Estime_LB2
SEXP C_Estime_LB2(SEXP s, int istart, int iend, double WS, double NEIGH, double sig);
RcppExport SEXP _Rnmr1D_C_Estime_LB2(SEXP sSEXP, SEXP istartSEXP, SEXP iendSEXP, SEXP WSSEXP, SEXP NEIGHSEXP, SEXP sigSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type s(sSEXP);
    Rcpp::traits::input_parameter< int >::type istart(istartSEXP);
    Rcpp::traits::input_parameter< int >::type iend(iendSEXP);
    Rcpp::traits::input_parameter< double >::type WS(WSSEXP);
    Rcpp::traits::input_parameter< double >::type NEIGH(NEIGHSEXP);
    Rcpp::traits::input_parameter< double >::type sig(sigSEXP);
    rcpp_result_gen = Rcpp::wrap(C_Estime_LB2(s, istart, iend, WS, NEIGH, sig));
    return rcpp_result_gen;
END_RCPP
}
// C_noise_estimate
SEXP C_noise_estimate(SEXP x, int n1, int n2, int flg);
RcppExport SEXP _Rnmr1D_C_noise_estimate(SEXP xSEXP, SEXP n1SEXP, SEXP n2SEXP, SEXP flgSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< int >::type n2(n2SEXP);
    Rcpp::traits::input_parameter< int >::type flg(flgSEXP);
    rcpp_result_gen = Rcpp::wrap(C_noise_estimate(x, n1, n2, flg));
    return rcpp_result_gen;
END_RCPP
}
// C_spec_ref_interval
SEXP C_spec_ref_interval(SEXP x, int istart, int iend, IntegerVector v);
RcppExport SEXP _Rnmr1D_C_spec_ref_interval(SEXP xSEXP, SEXP istartSEXP, SEXP iendSEXP, SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type istart(istartSEXP);
    Rcpp::traits::input_parameter< int >::type iend(iendSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(C_spec_ref_interval(x, istart, iend, v));
    return rcpp_result_gen;
END_RCPP
}
// C_spec_ref
SEXP C_spec_ref(SEXP x, IntegerVector v);
RcppExport SEXP _Rnmr1D_C_spec_ref(SEXP xSEXP, SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(C_spec_ref(x, v));
    return rcpp_result_gen;
END_RCPP
}
// C_MedianSpec
SEXP C_MedianSpec(SEXP x);
RcppExport SEXP _Rnmr1D_C_MedianSpec(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(C_MedianSpec(x));
    return rcpp_result_gen;
END_RCPP
}
// C_Derive1
SEXP C_Derive1(SEXP v);
RcppExport SEXP _Rnmr1D_C_Derive1(SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(C_Derive1(v));
    return rcpp_result_gen;
END_RCPP
}
// C_Derive
SEXP C_Derive(SEXP x);
RcppExport SEXP _Rnmr1D_C_Derive(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(C_Derive(x));
    return rcpp_result_gen;
END_RCPP
}
// C_Integre
SEXP C_Integre(SEXP x, int istart, int iend);
RcppExport SEXP _Rnmr1D_C_Integre(SEXP xSEXP, SEXP istartSEXP, SEXP iendSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type istart(istartSEXP);
    Rcpp::traits::input_parameter< int >::type iend(iendSEXP);
    rcpp_result_gen = Rcpp::wrap(C_Integre(x, istart, iend));
    return rcpp_result_gen;
END_RCPP
}
// C_segment_shifts
SEXP C_segment_shifts(SEXP x, int idx_vref, int decal_max, int istart, int iend, IntegerVector v);
RcppExport SEXP _Rnmr1D_C_segment_shifts(SEXP xSEXP, SEXP idx_vrefSEXP, SEXP decal_maxSEXP, SEXP istartSEXP, SEXP iendSEXP, SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type idx_vref(idx_vrefSEXP);
    Rcpp::traits::input_parameter< int >::type decal_max(decal_maxSEXP);
    Rcpp::traits::input_parameter< int >::type istart(istartSEXP);
    Rcpp::traits::input_parameter< int >::type iend(iendSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(C_segment_shifts(x, idx_vref, decal_max, istart, iend, v));
    return rcpp_result_gen;
END_RCPP
}
// C_align_segment
int C_align_segment(SEXP x, SEXP s, int istart, int iend, int apodize, IntegerVector v);
RcppExport SEXP _Rnmr1D_C_align_segment(SEXP xSEXP, SEXP sSEXP, SEXP istartSEXP, SEXP iendSEXP, SEXP apodizeSEXP, SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type x(xSEXP);
    Rcpp::traits::input_parameter< SEXP >::type s(sSEXP);
    Rcpp::traits::input_parameter< int >::type istart(istartSEXP);
    Rcpp::traits::input_parameter< int >::type iend(iendSEXP);
    Rcpp::traits::input_parameter< int >::type apodize(apodizeSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(C_align_segment(x, s, istart, iend, apodize, v));
    return rcpp_result_gen;
END_RCPP
}
// C_noise_estimation
double C_noise_estimation(SEXP x, int n1, int n2);
RcppExport SEXP _Rnmr1D_C_noise_estimation(SEXP xSEXP, SEXP n1SEXP, SEXP n2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< int >::type n2(n2SEXP);
    rcpp_result_gen = Rcpp::wrap(C_noise_estimation(x, n1, n2));
    return rcpp_result_gen;
END_RCPP
}
// C_aibin_buckets
SEXP C_aibin_buckets(SEXP x, SEXP b, SEXP v, SEXP l, int n1, int n2);
RcppExport SEXP _Rnmr1D_C_aibin_buckets(SEXP xSEXP, SEXP bSEXP, SEXP vSEXP, SEXP lSEXP, SEXP n1SEXP, SEXP n2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type x(xSEXP);
    Rcpp::traits::input_parameter< SEXP >::type b(bSEXP);
    Rcpp::traits::input_parameter< SEXP >::type v(vSEXP);
    Rcpp::traits::input_parameter< SEXP >::type l(lSEXP);
    Rcpp::traits::input_parameter< int >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< int >::type n2(n2SEXP);
    rcpp_result_gen = Rcpp::wrap(C_aibin_buckets(x, b, v, l, n1, n2));
    return rcpp_result_gen;
END_RCPP
}
// C_SDL_convolution
SEXP C_SDL_convolution(SEXP x, SEXP y, double sigma);
RcppExport SEXP _Rnmr1D_C_SDL_convolution(SEXP xSEXP, SEXP ySEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type x(xSEXP);
    Rcpp::traits::input_parameter< SEXP >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(C_SDL_convolution(x, y, sigma));
    return rcpp_result_gen;
END_RCPP
}
// C_erva_buckets
SEXP C_erva_buckets(SEXP x, SEXP b, SEXP v, SEXP l, int n1, int n2);
RcppExport SEXP _Rnmr1D_C_erva_buckets(SEXP xSEXP, SEXP bSEXP, SEXP vSEXP, SEXP lSEXP, SEXP n1SEXP, SEXP n2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type x(xSEXP);
    Rcpp::traits::input_parameter< SEXP >::type b(bSEXP);
    Rcpp::traits::input_parameter< SEXP >::type v(vSEXP);
    Rcpp::traits::input_parameter< SEXP >::type l(lSEXP);
    Rcpp::traits::input_parameter< int >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< int >::type n2(n2SEXP);
    rcpp_result_gen = Rcpp::wrap(C_erva_buckets(x, b, v, l, n1, n2));
    return rcpp_result_gen;
END_RCPP
}
// C_spectra_integrate
SEXP C_spectra_integrate(SEXP x, int istart, int iend);
RcppExport SEXP _Rnmr1D_C_spectra_integrate(SEXP xSEXP, SEXP istartSEXP, SEXP iendSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type istart(istartSEXP);
    Rcpp::traits::input_parameter< int >::type iend(iendSEXP);
    rcpp_result_gen = Rcpp::wrap(C_spectra_integrate(x, istart, iend));
    return rcpp_result_gen;
END_RCPP
}
// C_buckets_integrate
SEXP C_buckets_integrate(SEXP x, SEXP b, int mode);
RcppExport SEXP _Rnmr1D_C_buckets_integrate(SEXP xSEXP, SEXP bSEXP, SEXP modeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type x(xSEXP);
    Rcpp::traits::input_parameter< SEXP >::type b(bSEXP);
    Rcpp::traits::input_parameter< int >::type mode(modeSEXP);
    rcpp_result_gen = Rcpp::wrap(C_buckets_integrate(x, b, mode));
    return rcpp_result_gen;
END_RCPP
}
// C_all_buckets_integrate
SEXP C_all_buckets_integrate(SEXP x, SEXP b, int mode);
RcppExport SEXP _Rnmr1D_C_all_buckets_integrate(SEXP xSEXP, SEXP bSEXP, SEXP modeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type x(xSEXP);
    Rcpp::traits::input_parameter< SEXP >::type b(bSEXP);
    Rcpp::traits::input_parameter< int >::type mode(modeSEXP);
    rcpp_result_gen = Rcpp::wrap(C_all_buckets_integrate(x, b, mode));
    return rcpp_result_gen;
END_RCPP
}
// C_maxval_buckets
SEXP C_maxval_buckets(SEXP x, SEXP b);
RcppExport SEXP _Rnmr1D_C_maxval_buckets(SEXP xSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type x(xSEXP);
    Rcpp::traits::input_parameter< SEXP >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(C_maxval_buckets(x, b));
    return rcpp_result_gen;
END_RCPP
}
// C_ppmIntMax_buckets
SEXP C_ppmIntMax_buckets(SEXP x, SEXP b);
RcppExport SEXP _Rnmr1D_C_ppmIntMax_buckets(SEXP xSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type x(xSEXP);
    Rcpp::traits::input_parameter< SEXP >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(C_ppmIntMax_buckets(x, b));
    return rcpp_result_gen;
END_RCPP
}
// C_buckets_CSN_normalize
SEXP C_buckets_CSN_normalize(SEXP b);
RcppExport SEXP _Rnmr1D_C_buckets_CSN_normalize(SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(C_buckets_CSN_normalize(b));
    return rcpp_result_gen;
END_RCPP
}
// C_estime_sd
double C_estime_sd(SEXP x, int cut);
RcppExport SEXP _Rnmr1D_C_estime_sd(SEXP xSEXP, SEXP cutSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type cut(cutSEXP);
    rcpp_result_gen = Rcpp::wrap(C_estime_sd(x, cut));
    return rcpp_result_gen;
END_RCPP
}
// ajustBL
SEXP ajustBL(SEXP x, int flg);
RcppExport SEXP _Rnmr1D_ajustBL(SEXP xSEXP, SEXP flgSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type flg(flgSEXP);
    rcpp_result_gen = Rcpp::wrap(ajustBL(x, flg));
    return rcpp_result_gen;
END_RCPP
}
// C_corr_spec_re
SEXP C_corr_spec_re(SEXP l);
RcppExport SEXP _Rnmr1D_C_corr_spec_re(SEXP lSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type l(lSEXP);
    rcpp_result_gen = Rcpp::wrap(C_corr_spec_re(l));
    return rcpp_result_gen;
END_RCPP
}
// Fmin
double Fmin(SEXP par, SEXP re, SEXP im, int blphc, double B, int flg);
RcppExport SEXP _Rnmr1D_Fmin(SEXP parSEXP, SEXP reSEXP, SEXP imSEXP, SEXP blphcSEXP, SEXP BSEXP, SEXP flgSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type par(parSEXP);
    Rcpp::traits::input_parameter< SEXP >::type re(reSEXP);
    Rcpp::traits::input_parameter< SEXP >::type im(imSEXP);
    Rcpp::traits::input_parameter< int >::type blphc(blphcSEXP);
    Rcpp::traits::input_parameter< double >::type B(BSEXP);
    Rcpp::traits::input_parameter< int >::type flg(flgSEXP);
    rcpp_result_gen = Rcpp::wrap(Fmin(par, re, im, blphc, B, flg));
    return rcpp_result_gen;
END_RCPP
}
// Fentropy
double Fentropy(SEXP par, SEXP re, SEXP im, int blphc, int neigh, double B, double Gamma);
RcppExport SEXP _Rnmr1D_Fentropy(SEXP parSEXP, SEXP reSEXP, SEXP imSEXP, SEXP blphcSEXP, SEXP neighSEXP, SEXP BSEXP, SEXP GammaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type par(parSEXP);
    Rcpp::traits::input_parameter< SEXP >::type re(reSEXP);
    Rcpp::traits::input_parameter< SEXP >::type im(imSEXP);
    Rcpp::traits::input_parameter< int >::type blphc(blphcSEXP);
    Rcpp::traits::input_parameter< int >::type neigh(neighSEXP);
    Rcpp::traits::input_parameter< double >::type B(BSEXP);
    Rcpp::traits::input_parameter< double >::type Gamma(GammaSEXP);
    rcpp_result_gen = Rcpp::wrap(Fentropy(par, re, im, blphc, neigh, B, Gamma));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_Rnmr1D_lorentz", (DL_FUNC) &_Rnmr1D_lorentz, 4},
    {"_Rnmr1D_gauss", (DL_FUNC) &_Rnmr1D_gauss, 4},
    {"_Rnmr1D_pvoigt", (DL_FUNC) &_Rnmr1D_pvoigt, 5},
    {"_Rnmr1D_C_fSavGol", (DL_FUNC) &_Rnmr1D_C_fSavGol, 4},
    {"_Rnmr1D_C_FilterbyWT", (DL_FUNC) &_Rnmr1D_C_FilterbyWT, 4},
    {"_Rnmr1D_C_FilterbyThreshold", (DL_FUNC) &_Rnmr1D_C_FilterbyThreshold, 4},
    {"_Rnmr1D_C_Lorentz", (DL_FUNC) &_Rnmr1D_C_Lorentz, 5},
    {"_Rnmr1D_C_PVoigt", (DL_FUNC) &_Rnmr1D_C_PVoigt, 6},
    {"_Rnmr1D_C_OneVoigt", (DL_FUNC) &_Rnmr1D_C_OneVoigt, 3},
    {"_Rnmr1D_C_peakFinder", (DL_FUNC) &_Rnmr1D_C_peakFinder, 5},
    {"_Rnmr1D_C_peakOptimize", (DL_FUNC) &_Rnmr1D_C_peakOptimize, 4},
    {"_Rnmr1D_C_peakFiltering", (DL_FUNC) &_Rnmr1D_C_peakFiltering, 3},
    {"_Rnmr1D_C_specModel", (DL_FUNC) &_Rnmr1D_C_specModel, 3},
    {"_Rnmr1D_C_MyFuncTest2", (DL_FUNC) &_Rnmr1D_C_MyFuncTest2, 3},
    {"_Rnmr1D_SDL", (DL_FUNC) &_Rnmr1D_SDL, 2},
    {"_Rnmr1D_C_write_pack", (DL_FUNC) &_Rnmr1D_C_write_pack, 4},
    {"_Rnmr1D_C_read_pack", (DL_FUNC) &_Rnmr1D_C_read_pack, 1},
    {"_Rnmr1D_C_GlobSeg", (DL_FUNC) &_Rnmr1D_C_GlobSeg, 3},
    {"_Rnmr1D_lowpass1", (DL_FUNC) &_Rnmr1D_lowpass1, 2},
    {"_Rnmr1D_WinMoy", (DL_FUNC) &_Rnmr1D_WinMoy, 3},
    {"_Rnmr1D_Smooth", (DL_FUNC) &_Rnmr1D_Smooth, 2},
    {"_Rnmr1D_fitLines", (DL_FUNC) &_Rnmr1D_fitLines, 4},
    {"_Rnmr1D_C_Estime_LB", (DL_FUNC) &_Rnmr1D_C_Estime_LB, 6},
    {"_Rnmr1D_C_Estime_LB2", (DL_FUNC) &_Rnmr1D_C_Estime_LB2, 6},
    {"_Rnmr1D_C_noise_estimate", (DL_FUNC) &_Rnmr1D_C_noise_estimate, 4},
    {"_Rnmr1D_C_spec_ref_interval", (DL_FUNC) &_Rnmr1D_C_spec_ref_interval, 4},
    {"_Rnmr1D_C_spec_ref", (DL_FUNC) &_Rnmr1D_C_spec_ref, 2},
    {"_Rnmr1D_C_MedianSpec", (DL_FUNC) &_Rnmr1D_C_MedianSpec, 1},
    {"_Rnmr1D_C_Derive1", (DL_FUNC) &_Rnmr1D_C_Derive1, 1},
    {"_Rnmr1D_C_Derive", (DL_FUNC) &_Rnmr1D_C_Derive, 1},
    {"_Rnmr1D_C_Integre", (DL_FUNC) &_Rnmr1D_C_Integre, 3},
    {"_Rnmr1D_C_segment_shifts", (DL_FUNC) &_Rnmr1D_C_segment_shifts, 6},
    {"_Rnmr1D_C_align_segment", (DL_FUNC) &_Rnmr1D_C_align_segment, 6},
    {"_Rnmr1D_C_noise_estimation", (DL_FUNC) &_Rnmr1D_C_noise_estimation, 3},
    {"_Rnmr1D_C_aibin_buckets", (DL_FUNC) &_Rnmr1D_C_aibin_buckets, 6},
    {"_Rnmr1D_C_SDL_convolution", (DL_FUNC) &_Rnmr1D_C_SDL_convolution, 3},
    {"_Rnmr1D_C_erva_buckets", (DL_FUNC) &_Rnmr1D_C_erva_buckets, 6},
    {"_Rnmr1D_C_spectra_integrate", (DL_FUNC) &_Rnmr1D_C_spectra_integrate, 3},
    {"_Rnmr1D_C_buckets_integrate", (DL_FUNC) &_Rnmr1D_C_buckets_integrate, 3},
    {"_Rnmr1D_C_all_buckets_integrate", (DL_FUNC) &_Rnmr1D_C_all_buckets_integrate, 3},
    {"_Rnmr1D_C_maxval_buckets", (DL_FUNC) &_Rnmr1D_C_maxval_buckets, 2},
    {"_Rnmr1D_C_ppmIntMax_buckets", (DL_FUNC) &_Rnmr1D_C_ppmIntMax_buckets, 2},
    {"_Rnmr1D_C_buckets_CSN_normalize", (DL_FUNC) &_Rnmr1D_C_buckets_CSN_normalize, 1},
    {"_Rnmr1D_C_estime_sd", (DL_FUNC) &_Rnmr1D_C_estime_sd, 2},
    {"_Rnmr1D_ajustBL", (DL_FUNC) &_Rnmr1D_ajustBL, 2},
    {"_Rnmr1D_C_corr_spec_re", (DL_FUNC) &_Rnmr1D_C_corr_spec_re, 1},
    {"_Rnmr1D_Fmin", (DL_FUNC) &_Rnmr1D_Fmin, 6},
    {"_Rnmr1D_Fentropy", (DL_FUNC) &_Rnmr1D_Fentropy, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_Rnmr1D(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
