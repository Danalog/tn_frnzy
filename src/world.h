// tn_frnzy

/*
-----------------------------------------------------
Sound-Analysis Synthesis WORLD by M. Morise          |
                                                     |
In previous versions, fftw was not installed, but I  |
have install it for user simplicity                  |
-----------------------------------------------------
*/

// Currently known bugs
// decimateForF0 : Error just after the start and about 4 samples before the end.
//#include <fftsg.h>
#include "fftw/fftw3.h"
//#include <fft.h>

#include <stdlib.h>
#include <windows.h>
#include <math.h>

#define PI 3.1415926535897932384

// Only with windows
#pragma warning( disable : 4996 )

#pragma comment(lib, "libfftw3-3.lib")
#pragma comment(lib, "libfftw3f-3.lib")
#pragma comment(lib, "libfftw3l-3.lib")

#define MAX_FFT_LENGTH 2048
#define FLOOR_F0 90.0 // Prevent false F0 detection in low direction. optimistic that the original sound of UTAU is not so low. (Original value: 75.0)
#define DEFAULT_F0 500.0 // Size increased for cleaner consonants (Original value: 150.0)
#define LOW_LIMIT 65.0 // EFB-GT

// 71 is the lower limit that allows FFT length to be 2048 at fs: 44100.
// 70 Hz requires 4096 points.
// DEFAULT_F0 is a new feature. There is room for adjustment, but this is a tentative decision.

// FO Estimation with DIO : Distributed Inline filter Operation
void dio(double *x, int xLen, int fs, double framePeriod, 
		 double *timeAxis, double *f0);
int getSamplesForDIO(int fs, int xLen, double framePeriod);

// Spectral Envelope Estimation with STAR: Synchronous Technique and Adroit Restoration
int getFFTLengthForStar(int fs);

//void star(double *x, int xLen, int fs, double *timeAxis, double *f0,
//		  double **specgram);
//void getMinimumPhaseSpectrum(double *inputSpec, fftw_complex *spectrum, fftw_complex *cepstrum, int fftl);

// Aperiodicity index estimation method PLATINUM: Perfect Linear Aperiodicity Token Index Number Unified-Estimation Method
void pt100(double *x, int xLen, int fs, double *timeAxis, double *f0,  
		 double **residualSpecgram);

// Added from tn_fnds
int pt101(double *x, int xLen, int fs, double *timeAxis, double *f0, 
		 double ***residualSpecgram, int **residualSpecgramLength, int *residualSpecgramIndex);

// Added from tn_fnds
void PulseResidualWindow(double **residualSpecgram, int *residualSpecgramLength, int pCount);

// WORLD Synthesis
//void synthesis(double *f0, int tLen, double **specgram, double **residualSpecgram, int fftl, double framePeriod, int fs, 
//			   double *synthesisOut, int xLen);
void synthesisPt100(double *f0, int tLen, double **residualSpecgram, int fftl, double framePeriod, int fs, 
			   double *synthesisOut, int xLen);
//void getMinimumPhaseSpectrum(double *inputSpec, fftw_complex *spectrum, fftw_complex *cepstrum, int fftl);

// Added from tn_fnds
void synthesisPt101(double fixedDefault_f0, double *f0, int tLen, double **aperiodicity, int *ResidualSpecgramLength,
					int *fixedResidualSpecgramIndex, double *volume,
					int fftl, double framePeriod, int fs, double *synthesisOut, int xLen);
//------------------------------------------------------------------------------------
// Matlab Functions
double std2(double *x, int xLen);
void inv(double **r, int n, double **invr);
//void fftfilt(double *x, int xlen, double *h, int hlen, int fftl, double *y);
float randn(void);
void histc(double *x, int xLen, double *y, int yLen, int *index);
void interp1(double *t, double *y, int iLen, double *t1, int oLen, double *y1);
long decimateForF0(double *x, int xLen, double *y, int r);
void filterForDecimate(double *x, int xLen, double *y, int r);
int roundi(double x);
void diff(double *x, int xLength, double *ans);
void interp1Q(double x, double shift, double *y, int xLength, double *xi, int xiLength, double *ans);

