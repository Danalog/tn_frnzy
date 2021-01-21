
/*
Vocal Synthesis Engine tn_frnzy for UTAU, based on Zteer's tn_fnds

All comments were translated using DeepL
I decided on taking up the task of maintaining this resampler

you can find it at: https://github.com/Danalog/tn_frnzy

-- The following are comments from Zteer's tn_fnds

This program is based on Masamasa Morise's UTAU synthesis engine for WORLD, Eternal Force Blisampler. 
You can customize "Gentry Weeps: They Die, I Die" (EFB-GW) to create continuous sound and
consonant speed, and some flags.
I also used Ameya/Shobu's world4utau source code to create this.

 --- The following are comments from the original

Eternal Force Blissampler Gentry Weeps ~The other person dies, I die too~.
This is a UTAU synthesis engine for WORLD 0.0.4.
This program is not the same as WORLD.
The platinum number indicates the purity in thousandths, and platinum is recognized as platinum if it is above 850.
Therefore, Pt100 is something else (a new feature in this program) that cannot be considered platinum.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <windows.h>

#include "world.h"
#include "wavread.h"

#include <math.h>

// 13 Arguments:
// 1 - Input File (OK)
// 2 - Output File (OK)
// 3 - Scale (OK)
// 4 - Consonant Velocity
// 5 - Flag (OK)
// 6 - Offset
// 7 - Length Adjustment
// 8 - Consonant Part
// 9 - Blank
// 10 - Volume (OK)
// 11 - Modulation (OK)
// 12 - Tempo
// 13 - Pitch Bend

// Analysis shift amount [msec]
#define FRAMEPERIOD 2.0

//#pragma comment(lib, "winmm.lib")



void F0ToFile(double* f0, int tLen)
{
	FILE *file;
	int i;

	file = fopen("D:\\data\\YSS\\UTAU\\efb-gw20111228\\f0list.txt","w");
	for(i = 0; i < tLen; i++)
	{
		fprintf(file,"%d,%f\n",i ,f0[i]);
	}
	fclose(file);
}

void speqToFile(fftw_complex * spec, int fftl)
{
	FILE *file;
	int i;

	file = fopen("D:\\data\\YSS\\UTAU\\efb-gw20111228\\speqlist.txt","w");
	for(i = 0; i < fftl/2+1; i++)
	{
		fprintf(file,"%d,%f\n",i ,spec[i][0]);
	}
	fclose(file);
}


void createFinalPitch(double *f0, int tLen, double *pitchBend, int bLen, int signalLen, int offset, int fs, double tempo)
{
	int i;
	double *time1, *time2, *pitch;
	int pStep;
	pStep = (int)(60.0 / 96.0 / tempo * fs + 0.5);
	time1 = (double *)malloc(sizeof(double) * tLen);
	time2 = (double *)malloc(sizeof(double) * bLen);
	pitch = (double *)malloc(sizeof(double) * tLen);



	for(i = 0;i < tLen;i++) time1[i] = (double)i * FRAMEPERIOD;
    for(i = 0;i < bLen;i++) time2[i] = (double)i * pStep / (double)fs * 1000.0 + offset/1000.0;
	time2[0] = 0;
	interp1(time2, pitchBend, bLen, time1, tLen, pitch);

	for(i = (int)(offset*FRAMEPERIOD/1000);i < tLen;i++) f0[i] *= pitch[i];

//	for(i = 0;i < tLen;i+=10)
//	{
//		printf("%f\n", pitch[i]);
//	}

	free(time1); free(time2); free(pitch);
}

int base64decoderForUtau(char x, char y)
{
	int ans1, ans2, ans;

	if(x=='+') ans1 = 62;
	if(x=='/') ans1 = 63;
	if(x>='0' && x <= '9') ans1 = x+4;
	if(x>='A' && x <= 'Z') ans1 = x-65;
	if(x>='a' && x <= 'z') ans1 = x-71;

	if(y=='+') ans2 = 62;
	if(y=='/') ans2 = 63;
	if(y>='0' && y <= '9') ans2 = y+4;
	if(y>='A' && y <= 'Z') ans2 = y-65;
	if(y>='a' && y <= 'z') ans2 = y-71;

	ans = (ans1<<6) | ans2;
	if(ans >= 2048) ans -= 4096;
	return ans;
}

int getF0Contour(char *input, double *output)
{
	int i, j, count, length;
	i = 0;
	count = 0;
	double tmp;

	tmp = 0.0;
	while(input[i] != '\0')
	{
		if(input[i] == '#')
		{
			length = 0;
			for(j = i+1;input[j]!='#';j++)
			{
				length = length*10 + input[j]-'0';
			}
			i = j+1;
			for(j = 0;j < length;j++)
			{
				output[count++] = tmp;
			}
		}
		else
		{
			tmp = pow(2.0, (double)base64decoderForUtau(input[i], input[i+1]) / 1200.0);
			output[count++] = tmp;
			i+=2;
		}
	}

	return count;
}

// Ported from Ameya's world4utau.cpp
double getFreqAvg(double f0[], int tLen)
{
	int i, j;
	double value = 0, r;
	double p[6], q;
	double freq_avg = 0;
	double base_value = 0;
	for (i = 0; i < tLen; i++)
	{
		value = f0[i];
		if (value < 1000.0 && value > 55.0)
		{
			r = 1.0;
			// Heavier weighting for consecutive close values
			for (j = 0; j <= 5; j++)
			{
				if (i > j) {
					q = f0[i - j - 1] - value;
					p[j] = value / (value + q * q);
				} else {
					p[j] = 1/(1 + value);
				}
				r *= p[j];
			}
			freq_avg += value * r;
			base_value += r;
		}
	}
	if (base_value > 0) freq_avg /= base_value;
	return freq_avg;
}
// Ported from Ameya's world4utau.cpp
int get64(int c)
{
    if (c >= '0' && c <='9')
    {
        return c - '0' + 52;
    }
    else if (c >= 'A' && c <='Z')
    {
        return c - 'A';
    }
    else if (c >= 'a' && c <='z')
    {
        return c - 'a' + 26;
    }
    else if (c == '+')
    {
        return 62;
    }
    else if (c == '/')
    {
        return 63;
    }
    else
    {
        return 0;
    }
}
// Ported from Ameya's world4utau.cpp
int decpit(char *str, int *dst, int cnt)
{
	int len = 0;
	int i, n = 0;
	int k = 0, num, ii;
	if (str != NULL)
	{
		len = strlen(str);
		for (i = 0; i < len; i += 2)
		{
			if (str[i] == '#')
			{
				i++;
				sscanf(str + i, "%d", &num);
				for (ii = 0; ii < num && k < cnt; ii++) {
					dst[k++] = n;
				}
				while (str[i] != '#' && str[i] != 0) i++;
				i--;
			} 
			else
			{
				n = get64(str[i]) * 64 + get64(str[i + 1]);
				if (n > 2047) n -= 4096;
				if (k < cnt) {
					dst[k++] = n;
				}
			}
		}
	}
	return len;
}


void equalizingPicth(double *f0, int tLen, char *scaleParam, int modulationParam, int flag_t)
{
	int i;
	// Examine mean value
	double averageF0;
	double modulation;

	modulation = (double)modulationParam / 100.0;

	averageF0 = getFreqAvg(f0, tLen);

	int scale;
	int octave;
	double targetF0;
	int bias = 0;

	// Identify target scale
	if(scaleParam[1] == '#') bias = 1;

	switch(scaleParam[0])
	{
	case 'C':
		scale = -9+bias;
		break;
	case 'D':
		scale = -7+bias;
		break;
	case 'E':
		scale = -5;
		break;
	case 'F':
		scale = -4+bias;
		break;
	case 'G':
		scale = -2+bias;
		break;
	case 'A':
		scale = bias;
		break;
	case 'B':
		scale = 2;
		break;
	}
	octave = scaleParam[1+bias]-'0' - 4;
	targetF0 = 440 * pow(2.0,(double)octave) * pow(2.0, (double)scale/12.0);
	targetF0 *= pow(2, (double)flag_t/120);

	double tmp;
	
	if(averageF0 != 0.0)
	{
		for(i = 0;i < tLen;i++)
		{
			if(f0[i] != 0.0)
			{
				tmp = (f0[i]-averageF0)*modulation + averageF0;
				f0[i] = tmp * targetF0 / averageF0;
			}
		}
	}
	else
	{
		for(i = 0;i < tLen;i++)
		{
			if(f0[i] != 0.0)
			{
				f0[i] = targetF0;
			}
		}
	}
}

int stretchTime(double *f0, int tLen, int fftl, int *residualSpecgramIndex,
				 double *f02, int tLen2, int *residualSpecgramIndex2, int os, int st, int ed, int Length2, double vRatio, int mode)
{
	int i, k;
	int st2, ed2;

	st2 = min(tLen2, (int)((st-os) * vRatio + 0.5));  // Consonant frame after expansion and contraction
    ed2 = min(tLen2, (int)(Length2 + 0.5));     // Number of samples after synthesis
	// First half 
	for(i = 0;i < st2;i++)
	{
		k = max(0, min(tLen-1, int(i/vRatio) + os));
		f02[i] = f0[k];
		residualSpecgramIndex2[i] = residualSpecgramIndex[k];
	}
	// Second half (looped stretching)
	if(mode == 0)
	{
		i = st2;
		while(i < ed2)
		{
			for(k = st; k < ed - 2; k++)
			{
				if(i > ed2-1) break;
				f02[i] = f0[k];
				residualSpecgramIndex2[i] = residualSpecgramIndex[k];
				i++;
			}
			for(k = ed -1; k > st; k--)
				{
				if(i > ed2-1) break;
				f02[i] = f0[k];
				residualSpecgramIndex2[i] = residualSpecgramIndex[k];
				i++;
 			}
		}
	}
	else
	{
		// Second half (UTAU-Style stretching)
		if(ed2-st2 > ed-st) // Expansion
		{
			double ratio;
			ratio = (double)(ed-st)/(ed2-st2);
			for(i = st2;i < ed2; i++)
			{
				k = max(0, min(tLen-1, (int)((i - st2) * ratio + 0.5 + st)));
				f02[i] = f0[k];
				residualSpecgramIndex2[i] = residualSpecgramIndex[k];
			}
		}
		else
		{
			for(i = st2;i < ed2; i++)
			{
				k = st + (i - st2);
				f02[i] = f0[k];
				residualSpecgramIndex2[i] = residualSpecgramIndex[k];
			}
		}
	}

	return ed2;
}

void f0Lpf(double *f0, int tLen, int flag_d)
{
	int i;
	int addcount = 0;
	double addvalue = 0;
	double* newf0;
	newf0 = (double*)malloc(sizeof(double) * tLen);
		for(i = 0; i < min(tLen-1, flag_d); i++) 
	{
		if(f0[i] != 0.0)
		{
			addvalue += f0[i];
			addcount += 1;
		}
	}
	for(i = 0; i < tLen; i++)
	{	
			if(i - flag_d -1 >= 0)
		{
			if(f0[i - flag_d -1] != 0.0)
			{
				addvalue -= f0[i - flag_d -1];
				addcount -= 1;
			}
		}
		if(i + flag_d <= tLen - 1)
		{
			if(f0[i + flag_d] != 0.0)
			{
				addvalue += f0[i + flag_d];
				addcount += 1;
			}
		}
		if(f0[i] != 0)
		{
			newf0[i] = addvalue / addcount; 
		}
		else
		{
			newf0[i] = 0.0;
		}
	}
	for(i = 0; i < tLen; i++) f0[i] = newf0[i];
}

// Create spectrum for equalization
void createWaveSpec(double *x, int xLen, int fftl, int equLen, fftw_complex **waveSpecgram)
{
	int i, j;

	double *waveBuff;
	fftw_plan			wave_f_fft;				// Set FFT
	fftw_complex		*waveSpec;	// Spectral
	waveBuff = (double *)malloc(sizeof(double) * fftl);
	waveSpec = (fftw_complex *)malloc(sizeof(fftw_complex) * fftl);
	wave_f_fft = fftw_plan_dft_r2c_1d(fftl, waveBuff, waveSpec, FFTW_ESTIMATE);	

	int offset;

	for(i = 0;i < equLen;i++)
	{
		offset = i * fftl / 2;
		// Copy Data
		for(j = 0;j < fftl; j++) waveBuff[j] = x[offset + j] * 
										(0.5 - 0.5 * cos(2.0*PI*(double)j/(double)fftl));// Hang window

		// Execute FFT
		fftw_execute(wave_f_fft);

		// Store spectrogram
		for(j = 0;j < fftl/2+1; j++)
		{
			waveSpecgram[i][j][0] = waveSpec[j][0];
			waveSpecgram[i][j][1] = waveSpec[j][1];
		}
	}

	fftw_destroy_plan(wave_f_fft);
	free(waveBuff);
	free(waveSpec);

}

// Reconstruct waveforms from spectrogram
void rebuildWave(double *x, int xLen, int fftl, int equLen, fftw_complex **waveSpecgram)
{
	int i, j;
	double *waveBuff;
	fftw_plan			wave_i_fft;				// Set FFT
	fftw_complex		*waveSpec;	// Spectral
	waveBuff = (double *)malloc(sizeof(double) * fftl);
	waveSpec = (fftw_complex *)malloc(sizeof(fftw_complex) * fftl);
	wave_i_fft = fftw_plan_dft_c2r_1d(fftl, waveSpec, waveBuff, FFTW_ESTIMATE);	

	int offset;
	for(i = 0;i < xLen;i++) x[i] = 0;

	for(i = 0;i < equLen;i++)
	{
		offset = i * fftl / 2;

		// Store spectrogram
		for(j = 0;j < fftl/2+1; j++)
		{
			waveSpec[j][0] = waveSpecgram[i][j][0];
			waveSpec[j][1] = waveSpecgram[i][j][1];
		}


		// Execute FFT
		fftw_execute(wave_i_fft);

		for(j = 0;j < fftl; j++) waveBuff[j] /= fftl;

		// Copy Data
		for(j = 0;j < fftl; j++) x[offset + j]  += waveBuff[j]; 

	}

	fftw_destroy_plan(wave_i_fft);
	free(waveBuff);
	free(waveSpec);

}

/*
---------------
	FLAGS
---------------
*/

// B Flag (Breath)
void breath2(double *f0, int tLen, int fs, double *x, int xLen, fftw_complex **waveSpecgram,int equLen, int fftl, int flag_B)
{
	int i, j;

	// Preparation of FFT noise
	double *noiseData;
	double *noiseBuff;
	double *noise;
	fftw_plan			noise_f_fft;				// Set FFT
	fftw_plan			noise_i_fft;				// Set FFT
	fftw_complex		*noiseSpec;	// Spectral

	noiseData = (double *)malloc(sizeof(double) * xLen);
	//for(i=0;i < xLen; i++) noiseData[i] = (double)rand()/(RAND_MAX+1) - 0.5;
	for(i=0;i < xLen; i++) noiseData[i] = (double)rand()/((double)RAND_MAX+1) - 0.5;
	noise = (double *)malloc(sizeof(double) * xLen);
	for(i=0;i < xLen; i++) noise[i] = 0.0;
//	for(i=0;i < xLen; i++) noiseData[i] *= noiseData[i] * (noiseData[i] < 0)? -1 : 1; // Tweak noise distribution
	noiseBuff = (double *)malloc(sizeof(double) * fftl);
	noiseSpec = (fftw_complex *)malloc(sizeof(fftw_complex) * fftl);
	noise_f_fft = fftw_plan_dft_r2c_1d(fftl, noiseBuff, noiseSpec, FFTW_ESTIMATE);	
	noise_i_fft = fftw_plan_dft_c2r_1d(fftl, noiseSpec, noiseBuff, FFTW_ESTIMATE);	

	// Prepare wavefft
	fftw_complex		*waveSpec;	// Spectral
	waveSpec = (fftw_complex *)malloc(sizeof(fftw_complex) * fftl);

	int offset;
	double volume;

	int SFreq, MFreq, EFreq;

	SFreq = (int)(fftl * 1500 / fs); // Breath start frequency
	MFreq = (int)(fftl * 5000 / fs); // Breath middle frequency
	EFreq = (int)(fftl * 20000 / fs); // BReath end frequency

	double nowIndex;
	int sIndex, eIndex;
	double nowF0;
	int specs, spece;
	double hs, he;
	int baion;

	for(i = 0; i < equLen; i++)
	{
		offset = i * fftl / 2;
		// Copy Data
		for(j = 0;j < fftl; j++) noiseBuff[j] = noiseData[offset + j] *
										(0.5 - 0.5*cos(2.0*PI*(double)j/(double)fftl)); // Hang window

		// Execute FFT
		fftw_execute(noise_f_fft);

		// Spectral envelope (difficult)
		for(j = 0;j < fftl/2+1; j++) waveSpec[j][0] = sqrt(waveSpecgram[i][j][0] * waveSpecgram[i][j][0] + waveSpecgram[i][j][1] * waveSpecgram[i][j][1]);
		for(j = 0;j < fftl/2+1; j++) waveSpec[j][0] = log10(waveSpec[j][0]+0.00000001);//対数化
		for(j = 0;j < fftl/2+1; j++) waveSpec[j][1] = waveSpec[j][0];

		nowIndex = max(0, min(tLen-1, (double)(offset + fftl / 2) / fs * 1000 / FRAMEPERIOD));
		sIndex = min(tLen -2, (int)nowIndex);
		eIndex = sIndex + 1;
		
		nowF0 = (f0[sIndex] == 0 && f0[eIndex] == 0) ?  DEFAULT_F0 :
				(f0[sIndex] == 0) ? f0[eIndex] :
				(f0[eIndex] == 0) ? f0[sIndex] : 
									(f0[eIndex] - f0[sIndex]) * (nowIndex - sIndex) + f0[sIndex];

		specs = 0;
		hs = 0.0;
		j = 0;
		baion = 1;
		spece = 0;
		for(baion = 1;spece != fftl/2+1;baion++)
		{
			spece = min(fftl/2+1, (int)((double)fftl / fs * nowF0 * baion + 0.5));
			he = waveSpec[spece][1];
			for(j = specs;j < spece;j++)
			{
				waveSpec[j][0] = (he-hs)/(spece-specs)*(j-specs)+hs;
			}
			specs = spece;
			hs = he;
		}

		for(j = 0;j < fftl/2+1; j++) waveSpec[j][0] = pow(10, waveSpec[j][0]); // Amplification

		// Transform noise spectrogram
		for(j = 0;j < SFreq; j++)
		{
			noiseSpec[j][0] = 0.0;
			noiseSpec[j][1] = 0.0;
		}

		for(;j < MFreq; j++)
		{
			volume = waveSpec[j][0] * (0.5 - 0.5 * cos(PI * (j - SFreq) / (double)(MFreq - SFreq)));
			noiseSpec[j][0] *= volume;
			noiseSpec[j][1] *= volume;
		}
		for(;j < EFreq; j++)
		{
			volume = waveSpec[j][0] * (0.5 - 0.5 * cos(PI + PI * (j - MFreq) / (double)(EFreq - MFreq)));
			noiseSpec[j][0] *= volume;
			noiseSpec[j][1] *= volume;
		}

		for(;j < fftl/2+1; j++)
		{
			noiseSpec[j][0] = 0.0;
			noiseSpec[j][1] = 0.0;
		}
		
		noiseSpec[0][1] = 0.0;
		noiseSpec[fftl/2][1] = 0.0;

		// Invert FFT
		fftw_execute(noise_i_fft);
		for(j = 0;j < fftl; j++) noiseBuff[j] /= fftl;
		
		// Hang window
	//	for(j = 0;j < fftl; j++) noiseBuff[j] *= 0.5 - 0.5*cos(2.0*PI*(double)j/(double)fftl);

		// Add Noise
		// TODO: try messing with this
		for(j = 0;j < fftl; j++)
		{
			noise[offset + j] += noiseBuff[j] * 0.2;
		}
	}
	
	// Synthesize noise
	double noiseRatio = max(0, (double)(flag_B - 50) / 50.0);
	double waveRatio = 1 - noiseRatio;
	for(i = 0;i < xLen;i++) x[i] = x[i] * waveRatio + noise[i] * noiseRatio;

	// Post-Processing
	fftw_destroy_plan(noise_f_fft);
	fftw_destroy_plan(noise_i_fft);
	free(noise);
	free(noiseData);
	free(noiseBuff);
	free(noiseSpec);
	free(waveSpec);
}

// O Flag (Brightness)
void Opening(double *f0, int tLen, int fs, fftw_complex **waveSpecgram,int equLen, int fftl, int flag_O)
{
	int i, j;
	double opn = (double) flag_O / 100.0;
	int sFreq = (int)(fftl * 500 / fs); // Control frequency 1
	int eFreq = (int)(fftl * 2000 / fs); // Control frequency 2
	double sRatio = -10.0; // Amplitude multiplier of control frequency 1 (dB)
	double eRatio = 10.0; // Amplitude multiplier of control frequency 2 (dB)

	// Create volume map for each frequency
	double volume;
	double *volumeMap;
	volumeMap = (double *)malloc(sizeof(double) * fftl/2+1);

	volume = pow(10, sRatio * opn / 20);
	for(j = 0;j < sFreq;j++)
	{	
		volumeMap[j] = volume;
	}
	for(;j < eFreq;j++)
	{
		volume = pow(10, ((0.5+0.5*cos(PI+PI/(eFreq-sFreq)*(j-sFreq)))*(eRatio-sRatio)+sRatio) * opn / 20);
		volumeMap[j] = volume;
	}
	volume = pow(10, eRatio * opn / 20);
	for(;j < fftl/2+1;j++)
	{
		volumeMap[j] = volume;
	}

	// Change volume for each frequency
	int f0Frame;
	for(i = 0;i < equLen;i++)
	{
		f0Frame = max(0, min(tLen-1, (int)((double)((i+1) * fftl / 2) / fs * 1000 / FRAMEPERIOD + 0.5)));
		if(f0[f0Frame] == 0.0) continue;
		for(j = 0;j < fftl/2+1;j++)
		{	
			waveSpecgram[i][j][0] *= volumeMap[j];
			waveSpecgram[i][j][1] *= volumeMap[j];
		}
	}

	free(volumeMap);
}

// b Flag (Voiceless consonant stress)
void consonantAmp2(double *f0, double *volume, int tLen, int flag_b)
{
	int i;
	int frameLen = 5; //Number of frames to smooth (before/after)
	int addCount = 0;
	double addVolume = 0;
	double ratio = (double) flag_b / 20.0; // Magnification x5 when b=100

	for(i = 0;i < min(tLen, frameLen+1); i++)
	{
		addCount++;
		addVolume += (f0[i] == 0) ? ratio : 0.0;
	}
	for(i = 0;i < tLen-1; i++)
	{
		volume[i] *= (addCount != 0) ? addVolume / addCount + 1.0 : 1.0;

		if(i >= frameLen)
		{
			addCount--;
			addVolume -= (f0[i-frameLen] == 0) ? ratio : 0.0;
		}
		if(i <= tLen-1-frameLen-1)
		{
			addCount++;
			addVolume += (f0[i+frameLen+1] == 0) ? ratio : 0.0;
		}
	}
}

// g Flag (Gender Factor)
void gFactor(int pCount, int fftl, double **residualSpecgram, int *residualSpecgramLength, double gRatio)
{
	int i, j;
    double position;
	int sindex, eindex;
	int NewLength;

	for(i = 0; i < pCount-1; i++)
	{
		if(residualSpecgramLength[i] == 0.0) continue;

		NewLength = max(0, min(fftl-1, (int)(residualSpecgramLength[i] / gRatio + 0.5)));
		if (gRatio>1)
		{
			for(j = 0;j < NewLength;j++)
			{
				position = min((double)residualSpecgramLength[i]-1.0001, (double)(j * gRatio));
				sindex = (int)position;
				eindex = sindex + 1;
				residualSpecgram[i][j] = residualSpecgram[i][sindex] + 
					             (double)(residualSpecgram[i][eindex] - residualSpecgram[i][sindex]) * 
								 (double)(position - sindex);
			}
		}
		else
		{
			for(j = NewLength-1;j >= 0;j--)
			{
				position = min((double)residualSpecgramLength[i]-1.0001, (double)(j * gRatio));
				sindex = (int)position;
				eindex = sindex + 1;
				residualSpecgram[i][j] = residualSpecgram[i][sindex] + 
					             (double)(residualSpecgram[i][eindex] - residualSpecgram[i][sindex]) * 
								 (double)(position - sindex);
			}
		}
		residualSpecgramLength[i] = NewLength;
	}
}

// Correct the period of the noise part changed by the g flag (adjust F0 to the period)
// Execute just before synthesisPt101 because F0 of the noise part will not be 0.0
void f0FixG(double *f0, int tLen2, double gRatio)
{
	int i;
	for(i = 0;i < tLen2;i++)
	{
		if(f0[i] == 0.0)
		{
			f0[i] = DEFAULT_F0 * gRatio;
		}
	}
}

// Add noise to F0
void f0Noise(double *f0, int tLen, double f0Rand)
{
	int i, j;
	int Pit = 1;
	double sRand, eRand;
	double NowRand;

	eRand = 0;
	i = 0;
	while(i <= tLen-1)
	{
		sRand = eRand;
	//	eRand = (double)rand()/(RAND_MAX+1) * f0Rand * 2  - f0Rand;
	//	eRand = (double)rand()/(RAND_MAX+1) * -f0Rand;
		eRand = (double)(int)(rand()/(RAND_MAX/3)-1) * f0Rand;
		for(j = 0;(j < Pit) && (i+j <= tLen-1); j++)
		{
			if(f0[i+j] != 0.0)
			{
				NowRand = (eRand - sRand) / Pit * j + sRand;
				f0[i+j] *= pow(2, NowRand);
			}
		}
		i += j;
	}
}
// Convert frequency to pitch

double FrqToPit(double Frq)
{
	return log(Frq / 220) * 1.44269504088896 * 1200 + 5700;
}

// A Flag (Maps volume to pitch change)
void autoVolume(double *f0, int tLen, int fs, double *volume, int flag_A)
{
	int i;
	
	if(flag_A == 0)
	{
		for(i = 0;i < tLen; i++) volume[i] = 1.0;
		return;
	}

	double AutoPow;
	for(i = 0;i < tLen-1; i++)
	{
		if(f0[i] == 0.0)
		{
			volume[i] = 1.0;	
			continue;
		}

		if (f0[i+1] != 0.0)	
		{
			AutoPow = (FrqToPit(f0[i+1]) - FrqToPit(f0[i])) * (441 / (fs * FRAMEPERIOD)) * flag_A; 
			volume[i] = min(1.2, pow(2, AutoPow * 1));

			continue;
		}

		if(i > 0)
		{
			if(f0[i-1] != 0.0)
			{
				volume[i] = volume[i-1];
				continue;
			}
		}
		volume[i] = 1.0;
	}
	if(f0[tLen-1] != 0.0 && f0[tLen-2] != 0.0) volume[tLen-1] = volume[tLen-2];
}

int main(int argc, char *argv[])
{
	int i;

	double *x,*f0,*t,*y;
	double **residualSpecgram;
	int *residualSpecgramLength;
	int *residualSpecgramIndex;
	int pCount;
	int fftl;

	int signalLen;
	int tLen;

	if(argc < 3) 
	{
		printf("error: the number of parameters is not correct.\n");
		return 0;
	}

	/*
	printf("argc:%d\n", argc);
	for(i = 0;i < argc-1;i++)
		printf("%s\n", argv[i]);
	//*/

	// Read flags
	char *cp;
	int flag_B = 50; // B Flag (Breath)
	if(argc > 5 && (cp = strchr(argv[5],'B')) != 0)
	{
		sscanf(cp+1, "%d", &flag_B);
		flag_B = max(0, min(100, flag_B));
	}

	int flag_b = 0; // b Flag (Voiceless consonant stress)
	if(argc > 5 && (cp = strchr(argv[5],'b')) != 0)
	{
		sscanf(cp+1, "%d", &flag_b);
		flag_b = max(0, min(100, flag_b));
	}

	int flag_t = 0; // t Flag (Pitch shift (in cents))
	if(argc > 5 && (cp = strchr(argv[5],'t')) != 0)
	{
		sscanf(cp+1, "%d", &flag_t);
	}

	double flag_g = 0.0; // g Flag (Gender Factor)
	double gRatio;
	if(argc > 5 && (cp = strchr(argv[5],'g')) != 0)
	{
		sscanf(cp+1, "%lf", &flag_g);
		if (flag_g>100) flag_g = 100;
		if (flag_g<-100) flag_g= -100;
	}
	gRatio = pow(10, -flag_g / 200);

	double flag_W = 0.0; // W Flag（Forced Frequency Setting）F<0 silent F=0 disabled 50>=F<=1000 Set to specified frequency  
	double f0Rand = 0;
	if(argc > 5 && (cp = strchr(argv[5],'W')) != 0)
	{
		sscanf(cp+1, "%lf", &flag_W);
		if (flag_W > 1000) flag_W = 1000;
		if ((flag_W <    50) && (flag_W >    0)){f0Rand =  flag_W / 50; flag_W = 0;}
		if (flag_W <    0) flag_W = -1;
	}
	int flag_d = 5; // d Flag (Apply LPF to DIO F0 analysis results 0~20 def 5)
	if(argc > 5 && (cp = strchr(argv[5],'d')) != 0)
	{
		sscanf(cp+1, "%d", &flag_d);
		flag_d = max(0, min(20, flag_d));
	}

	int flag_A = 0; // A Flag　(Map volume to pitch change)
	if(argc > 5 && (cp = strchr(argv[5],'A')) != 0)
	{
		sscanf(cp+1, "%d", &flag_A);
		flag_A = max(0, min(100, flag_A));
	}

	int flag_O = 0;// O Flag (Brightness)
	if(argc > 5 && (cp = strchr(argv[5],'O')) != 0)
	{
		sscanf(cp+1, "%d", &flag_O);
		flag_O = max(-100, min(100, flag_O));
	}

	int flag_e = 0;// e Flag (Stretching Method, 1=Looping 2=UTAU-Style stretching)
	if(argc > 5 && (cp = strchr(argv[5],'e')) != 0)
	{
		flag_e = 1;
	}



	FILE *fp;

	int offset;
	int edLengthMsec;
	offset = atoi(argv[6]);
	edLengthMsec = atoi(argv[9]);
	int stp = offset;

	int fs, nbit;
	//x = wavread(argv[1], &fs, &nbit, &signalLen);
	x = wavread(argv[1], &fs, &nbit, &signalLen, &offset, &edLengthMsec);

	if(x == NULL)
	{
		printf("error: The specified file does not exist.\n");
		return 0;
	}

	printf("File information\n");
	printf("Sampling : %d Hz %d Bit\n", fs, nbit);
	printf("Length %d [sample]\n", signalLen);
	printf("Length %f [sec]\n", (double)signalLen/(double)fs);

	// Pre-Calculate amount of samples in F0．
	tLen = getSamplesForDIO(fs, signalLen, FRAMEPERIOD);
	f0 = (double *)malloc(sizeof(double)*tLen);
	t  = (double *)malloc(sizeof(double)*tLen);
	// F0 Estimation with DIO
	DWORD elapsedTime;
	
	if(flag_W == 0) // W Flag (Forced F0 Setting)
	{
		printf("\nAnalysis\n");
		elapsedTime = timeGetTime();

		dio(x, signalLen, fs, FRAMEPERIOD, t, f0);
		printf("DIO: %d [msec]\n", timeGetTime() - elapsedTime);


//		F0ToFile(f0, tLen);
		//F0のLPF  
		if (flag_d !=0)
		{
			f0Lpf(f0, tLen, flag_d);
		}
//		F0ToFile(f0, tLen);

	}
	else
	{
		for(i = 0;i < tLen;i++)
		{
			f0[i] = (flag_W == -1) ? 0.0 : flag_W;
			t[i] = (double)i * FRAMEPERIOD/1000.0;
		}
	}
	
	fftl = getFFTLengthForStar(fs);

	// Analysis of Aperiodic Indicators
	elapsedTime = timeGetTime();
	residualSpecgramIndex = (int *)malloc(sizeof(int) * tLen);

	pCount = pt101(x, signalLen, fs, t, f0, &residualSpecgram, &residualSpecgramLength, residualSpecgramIndex);
	printf("PLATINUM: %d [msec]\n", timeGetTime() - elapsedTime);

// Apply g flag
	if(flag_g != 0)
	{
		 gFactor(pCount, fftl, residualSpecgram, residualSpecgramLength, gRatio);
	}

	// Hang window
	PulseResidualWindow(residualSpecgram, residualSpecgramLength, pCount);

	// Expansion and Contraction of time length
	int lengthMsec, stLengthMsec, /*edLengthMsec,*/ inputLengthMsec;
	double velocity;
	double vRatio;

	inputLengthMsec = (int)(tLen*FRAMEPERIOD);     // Usable length of original sound
	lengthMsec = atoi(argv[7]);                    // Requested length
	stLengthMsec = atoi(argv[8]);                  // Consonant
	velocity = (double)atoi(argv[4]);              // Consonant Velocity
	vRatio = pow(2.0, (1.0 - (velocity / 100.0))); // Consonant Expansion Rate

	// Memory allocation for control parameters
	double *fixedF0;
	int *fixedResidualSpecgramIndex;
	double *fixedVolume;         // Volume per frame

	int tLen2;

    tLen2 = (int)(0.5+(double)(lengthMsec  )/FRAMEPERIOD);

	fixedF0					= (double *) malloc(sizeof(double)   * tLen2);
	fixedResidualSpecgramIndex	= (int *) malloc(sizeof(int) * tLen2);
	fixedVolume	= (double *) malloc(sizeof(double) * tLen2);

	// Allocate memory for final waveform
	int signalLen2;
	signalLen2 = (int)((lengthMsec       )/1000.0*(double)fs);
	y  = (double *)malloc(sizeof(double)*signalLen2);
	for(i = 0;i < signalLen2;i++) y[i] = 0.0;
//	printf("length:%d, %f\n",signalLen2, (double)signalLen2/(double)fs*1000.0);
//	printf("%d, %d, %d\n",lengthMsec, offset, fs);


	// F0 Operation before synthesis (argument)
	equalizingPicth(f0, tLen, argv[3], atoi(argv[11]), flag_t );

	// Time Stretching
	int os, st, ed;
	os = offset;
	st = stLengthMsec + offset;
	ed = inputLengthMsec - edLengthMsec;

	tLen2 = stretchTime(f0, tLen, fftl, residualSpecgramIndex, 
			fixedF0, tLen2, fixedResidualSpecgramIndex,
			os/(int)FRAMEPERIOD, st/(int)FRAMEPERIOD, min(ed/(int)FRAMEPERIOD, tLen-1),
			lengthMsec, vRatio, flag_e);

	// Apply pitch bend. Uses processing from world4utau
	int *pit = NULL;
	double tempo = 120;
	int pLen = tLen2;
	int pStep = 256;
	if (argc > 13) 	
	{
		cp = argv[12];
		sscanf(cp + 1, "%lf", &tempo);
		pStep = (int)(60.0 / 96.0 / tempo * fs + 0.5);
		pLen = signalLen2 / pStep + 1;
		pit = (int*)malloc((pLen+1) * sizeof(int) );
		memset(pit, 0, (pLen+1) * sizeof(int));
		decpit(argv[13], pit, pLen);
	}
	else
	{
		pit = (int*)malloc((pLen+1) * sizeof(int));
		memset(pit, 0, (pLen+1) * sizeof(int));
	}

	double tmo;
	double u;
	int m;
	for (i = 0; i < tLen2; i++)
  	{
		tmo = FRAMEPERIOD * i;
		u = tmo * 0.001 * fs / pStep;
		m = (int)floor(u);
		u -= m;
		if (m >= pLen) m = pLen - 1;
		fixedF0[i] *= pow(2, (pit[m] * (1.0 - u) + pit[m + 1] * u) / 1200.0);
	}
//	createFinalPitch(fixedF0, tLen2, pitchBend, bLen, signalLen2, offset, fs, tempo);

	// W flag pitch noise I've been plotting to de-voice, but it's not working.
//	if(f0Rand != 0.0)
//	{
//		f0Noise(fixedF0, tLen2, f0Rand);

//	}
	// Apply A Flag
	autoVolume(fixedF0, tLen2, fs, fixedVolume, flag_A);

	// Apply b Flag
	if(flag_b != 0)
	{
		consonantAmp2(fixedF0, fixedVolume, tLen2, flag_b);
	}

	// Reconcile the period of the silent part mismatched by the g flag with F0.
	double fixedDefault_f0 = DEFAULT_F0 * gRatio;

	// Synthesis
	printf("\n[tn_frnzy] Synthesis\n");
	elapsedTime = timeGetTime();
	synthesisPt101(fixedDefault_f0, fixedF0, tLen2, residualSpecgram, residualSpecgramLength, fixedResidualSpecgramIndex,
		fixedVolume, fftl, FRAMEPERIOD, fs, y, signalLen2);

	printf("WORLD: %d [msec]\n", timeGetTime() - elapsedTime);

	// Equalization
	int equfftL = 1024; // Equalize FFT Length
	int equLen = (signalLen2 / (equfftL/2)) - 1; // Number of repititions
	fftw_complex **waveSpecgram;  // Spectral
	waveSpecgram = (fftw_complex **)malloc(sizeof(fftw_complex *) * equLen);
	for(i = 0;i < equLen;i++) waveSpecgram[i] = (fftw_complex *)malloc(sizeof(fftw_complex) * (equfftL/2+1));

	// Create Spectrogram
	if(flag_B > 50 || flag_O != 0)
	{
		createWaveSpec(y, signalLen2, equfftL, equLen, waveSpecgram);
	}

	// Apply O Flag
	if(flag_O != 0)
	{
		Opening(fixedF0, tLen2, fs, waveSpecgram, equLen, equfftL, flag_O);
	}

	// Reflect equalization result in waveform
	if(flag_O != 0)
	{
		rebuildWave(y, signalLen2, equfftL, equLen, waveSpecgram);
	}

	// B Flag
	if(flag_B > 50)
	{
		 breath2(fixedF0, tLen2, fs, y, signalLen2, waveSpecgram, equLen, equfftL, flag_B);
	}

	// Set Offset
//	signalLen2 = (int)((lengthMsec)/1000.0*(double)fs);

	// Exporting a file (not related to its content)
	char header[44];
	short *output;
	double maxAmp;
	output = (short *)malloc(sizeof(short) * signalLen2);
 
	// Amplitude Normalization
	maxAmp = 0.0;
	double volume;
	volume = (double)atoi(argv[10]) / 100.0;
	for(i = 0;i < signalLen2;i++) maxAmp = maxAmp < fabs(y[i]) ? fabs(y[i]) : maxAmp;
	for(i = 0;i < signalLen2;i++) output[i] = (short)(32768.0*(y[i]*0.5 * volume/maxAmp));

	fp = fopen(argv[1], "rb");
	fread(header, sizeof(char), 44, fp);
	fclose(fp);

	*((short int*)(&header[22])) = 1;		 //channels	 	2 	# of Channels
	*((int*)(&header[24])) = fs;			 //samplerate 	4 	samples/sec
	*((int*)(&header[28])) = fs * nbit / 8;	 //bytepersec 	4 	bytes/sec
	*((short int*)(&header[32])) = nbit / 8; //blockalign 	2 	bytes/block
	*((short int*)(&header[34])) = nbit;	 //bitswidth 	2 	bits/sample

	header[36] = 'd'; header[37] = 'a'; header[38] = 't'; header[39] = 'a';

	fp = fopen(argv[2],"wb");
	fwrite(header, sizeof(char), 44, fp);
	fwrite(output, sizeof(short), signalLen2, fp);
	fseek(fp, 40, SEEK_SET);
	signalLen2*=2;
	fwrite(&signalLen2, sizeof(int), 1, fp);
	fclose(fp);
	free(output);

	free(pit);
	free(x); free(t); free(f0); free(fixedF0); free(y);
	for(i = 0;i < pCount;i++)
	{
		free(residualSpecgram[i]);
	}
	free(residualSpecgram); 
	free(fixedResidualSpecgramIndex);
	free(fixedVolume);
	free(residualSpecgramIndex); 
	free(residualSpecgramLength); 

	for(i = 0;i < equLen;i++) free(waveSpecgram[i]);
	free(waveSpecgram);

	printf("Complete!\n");

	return 0;
}
