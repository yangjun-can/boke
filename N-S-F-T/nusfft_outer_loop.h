#pragma once
#ifndef NUSFFT_OUTER_LOOP_H
#define NUSFFT_OUTER_LOOP_H

#include <complex>
#include <complex.h>
//#include <cmath>
//#include <complex.h>
//#include <vector>
//#include <map>
using namespace std;

using std::complex;
#define Vec(a, b) std::vector<__typeof(*(a))> ((a), (a)+(b))

constexpr auto OPTIMIZE_FFTW = 0;
//fft.h
// allow easy change to float or long double
//#define USE_FLOAT
#define USE_DOUBLE

#ifdef USE_FLOAT
//typedef float complex complex_t;
typedef complex<float> complex_t;
typedef float real_t;
#define cexp cexpf
#define exp expf
#endif

#ifdef USE_DOUBLE
//typedef double complex complex_t;
typedef complex<double> complex_t;
typedef double real_t;
#endif

//computerfourier.h

//#define  WITH_COMB 0 
//extern bool WITH_COMB;
//extern bool ALGORITHM1;
//extern bool VERBOSE;
//extern bool TIMING;

//fftw.h
#ifdef USE_FLOAT

#define fftw_plan_dft_1d fftwf_plan_dft_1d 
#define fftw_plan fftwf_plan
#define fftw_execute fftwf_execute
#define fftw_destroy_plan fftwf_destroy_plan
#define fftw_free fftwf_free

#endif


#endif