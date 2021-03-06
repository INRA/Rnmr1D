# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

lorentz <- function(x, x0, s) {
    .Call('_Rnmr1D_lorentz', PACKAGE = 'Rnmr1D', x, x0, s)
}

C_fSavGol <- function(s, m, nl, nr) {
    .Call('_Rnmr1D_C_fSavGol', PACKAGE = 'Rnmr1D', s, m, nl, nr)
}

C_FilterbyWT <- function(s, type, threshold = 0.5, verbose = 0L) {
    .Call('_Rnmr1D_C_FilterbyWT', PACKAGE = 'Rnmr1D', s, type, threshold, verbose)
}

C_FilterbyThreshold <- function(s, wavelet, threshold = 0L) {
    .Call('_Rnmr1D_C_FilterbyThreshold', PACKAGE = 'Rnmr1D', s, wavelet, threshold)
}

C_Lorentz <- function(ppm, amp, x0, sigma) {
    .Call('_Rnmr1D_C_Lorentz', PACKAGE = 'Rnmr1D', ppm, amp, x0, sigma)
}

C_OneLorentz <- function(X, Y, par) {
    .Call('_Rnmr1D_C_OneLorentz', PACKAGE = 'Rnmr1D', X, Y, par)
}

C_MyFuncTest <- function(spec, ppmrange, filt = NULL, peaks = NULL, verbose = 1L) {
    .Call('_Rnmr1D_C_MyFuncTest', PACKAGE = 'Rnmr1D', spec, ppmrange, filt, peaks, verbose)
}

C_MyFuncTest2 <- function(spec, n1, n2) {
    .Call('_Rnmr1D_C_MyFuncTest2', PACKAGE = 'Rnmr1D', spec, n1, n2)
}

SDL <- function(x, Sigma) {
    .Call('_Rnmr1D_SDL', PACKAGE = 'Rnmr1D', x, Sigma)
}

C_write_pack <- function(x, pmin, pmax, ff) {
    invisible(.Call('_Rnmr1D_C_write_pack', PACKAGE = 'Rnmr1D', x, pmin, pmax, ff))
}

C_read_pack <- function(ff) {
    .Call('_Rnmr1D_C_read_pack', PACKAGE = 'Rnmr1D', ff)
}

C_GlobSeg <- function(v, dN, sig) {
    .Call('_Rnmr1D_C_GlobSeg', PACKAGE = 'Rnmr1D', v, dN, sig)
}

lowpass1 <- function(x, alpha) {
    .Call('_Rnmr1D_lowpass1', PACKAGE = 'Rnmr1D', x, alpha)
}

WinMoy <- function(v, n1, n2) {
    .Call('_Rnmr1D_WinMoy', PACKAGE = 'Rnmr1D', v, n1, n2)
}

Smooth <- function(v, n) {
    .Call('_Rnmr1D_Smooth', PACKAGE = 'Rnmr1D', v, n)
}

fitLines <- function(s, b, n1, n2) {
    invisible(.Call('_Rnmr1D_fitLines', PACKAGE = 'Rnmr1D', s, b, n1, n2))
}

C_Estime_LB <- function(s, istart, iend, WS, NEIGH, sig) {
    .Call('_Rnmr1D_C_Estime_LB', PACKAGE = 'Rnmr1D', s, istart, iend, WS, NEIGH, sig)
}

C_Estime_LB2 <- function(s, istart, iend, WS, NEIGH, sig) {
    .Call('_Rnmr1D_C_Estime_LB2', PACKAGE = 'Rnmr1D', s, istart, iend, WS, NEIGH, sig)
}

C_noise_estimate <- function(x, n1, n2, flg) {
    .Call('_Rnmr1D_C_noise_estimate', PACKAGE = 'Rnmr1D', x, n1, n2, flg)
}

C_spec_ref_interval <- function(x, istart, iend, v) {
    .Call('_Rnmr1D_C_spec_ref_interval', PACKAGE = 'Rnmr1D', x, istart, iend, v)
}

C_spec_ref <- function(x, v) {
    .Call('_Rnmr1D_C_spec_ref', PACKAGE = 'Rnmr1D', x, v)
}

C_MedianSpec <- function(x) {
    .Call('_Rnmr1D_C_MedianSpec', PACKAGE = 'Rnmr1D', x)
}

C_Derive1 <- function(v) {
    .Call('_Rnmr1D_C_Derive1', PACKAGE = 'Rnmr1D', v)
}

C_Derive <- function(x) {
    .Call('_Rnmr1D_C_Derive', PACKAGE = 'Rnmr1D', x)
}

C_Integre <- function(x, istart, iend) {
    .Call('_Rnmr1D_C_Integre', PACKAGE = 'Rnmr1D', x, istart, iend)
}

C_segment_shifts <- function(x, idx_vref, decal_max, istart, iend, v) {
    .Call('_Rnmr1D_C_segment_shifts', PACKAGE = 'Rnmr1D', x, idx_vref, decal_max, istart, iend, v)
}

C_align_segment <- function(x, s, istart, iend, v) {
    .Call('_Rnmr1D_C_align_segment', PACKAGE = 'Rnmr1D', x, s, istart, iend, v)
}

C_noise_estimation <- function(x, n1, n2) {
    .Call('_Rnmr1D_C_noise_estimation', PACKAGE = 'Rnmr1D', x, n1, n2)
}

C_aibin_buckets <- function(x, b, v, l, n1, n2) {
    .Call('_Rnmr1D_C_aibin_buckets', PACKAGE = 'Rnmr1D', x, b, v, l, n1, n2)
}

C_spectra_integrate <- function(x, istart, iend) {
    .Call('_Rnmr1D_C_spectra_integrate', PACKAGE = 'Rnmr1D', x, istart, iend)
}

C_buckets_integrate <- function(x, b, mode) {
    .Call('_Rnmr1D_C_buckets_integrate', PACKAGE = 'Rnmr1D', x, b, mode)
}

C_all_buckets_integrate <- function(x, b, mode) {
    .Call('_Rnmr1D_C_all_buckets_integrate', PACKAGE = 'Rnmr1D', x, b, mode)
}

C_maxval_buckets <- function(x, b) {
    .Call('_Rnmr1D_C_maxval_buckets', PACKAGE = 'Rnmr1D', x, b)
}

C_ppmIntMax_buckets <- function(x, b) {
    .Call('_Rnmr1D_C_ppmIntMax_buckets', PACKAGE = 'Rnmr1D', x, b)
}

C_buckets_CSN_normalize <- function(b) {
    .Call('_Rnmr1D_C_buckets_CSN_normalize', PACKAGE = 'Rnmr1D', b)
}

C_estime_sd <- function(x, cut) {
    .Call('_Rnmr1D_C_estime_sd', PACKAGE = 'Rnmr1D', x, cut)
}

ajustBL <- function(x, flg) {
    .Call('_Rnmr1D_ajustBL', PACKAGE = 'Rnmr1D', x, flg)
}

C_corr_spec_re <- function(l) {
    .Call('_Rnmr1D_C_corr_spec_re', PACKAGE = 'Rnmr1D', l)
}

Fmin <- function(par, re, im, blphc, B, flg = 0L) {
    .Call('_Rnmr1D_Fmin', PACKAGE = 'Rnmr1D', par, re, im, blphc, B, flg)
}

Fentropy <- function(par, re, im, blphc, neigh, B, Gamma) {
    .Call('_Rnmr1D_Fentropy', PACKAGE = 'Rnmr1D', par, re, im, blphc, neigh, B, Gamma)
}

C_SDL_convolution <- function(x, y, sigma) {
    .Call('_Rnmr1D_C_SDL_convolution', PACKAGE = 'Rnmr1D', x, y, sigma)
}

