// wavread

#include <windows.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "wavread.h"

//#pragma warning(disable:4996)

/* Port wavread function */
double * wavread(char* filename, int *fs, int *Nbit, int *waveLength, int *offset, int *endbr)
{
	FILE *fp;
	char dataCheck[5]; // More
	unsigned char forIntNumber[4];
	double tmp, signBias, zeroLine;
	short int channel;
	int quantizationByte;
	double *waveForm;

	dataCheck[4] = '\0'; // For string collation, insert a terminating character at the end.
	fp = fopen(filename, "rb");
	if(NULL == fp) 
	{
	  printf("Loading file %s failed.\n", filename);
		return NULL;
	}
	// Header check
	fread(dataCheck, sizeof(char), 4, fp); // "RIFF"
	if(0 != strcmp(dataCheck,"RIFF"))
	{
		fclose(fp);
		printf("Incorrect RIFF header\n");
		return NULL;
	}
	fseek(fp, 4, SEEK_CUR); // Skip 4 bytes
	fread(dataCheck, sizeof(char), 4, fp); // "WAVE"
	if(0 != strcmp(dataCheck,"WAVE"))
	{
		fclose(fp);
		printf("Incorrect WAVE header\n");
		return NULL;
	}
	fread(dataCheck, sizeof(char), 4, fp); // "fmt "
	if(0 != strcmp(dataCheck,"fmt "))
	{
		fclose(fp);
		printf("Incorrect fmt header\n");
		return NULL;
	}
	fread(dataCheck, sizeof(char), 4, fp); //1 0 0 0
	if(!(16 == dataCheck[0] && 0 == dataCheck[1] && 0 == dataCheck[2] && 0 == dataCheck[3]))
	{
		fclose(fp);
		printf("Incorrect fmt header (2)\n");
		return NULL;
	}
	fread(dataCheck, sizeof(char), 2, fp); //1 0
	if(!(1 == dataCheck[0] && 0 == dataCheck[1]))
	{
		fclose(fp);
		printf("Incorrect format ID\n");
		return NULL;
	}
	/*
	fread(dataCheck, sizeof(char), 2, fp); //1 0
	if(!(1 == dataCheck[0] && 0 == dataCheck[1]))
	{
		fclose(fp);
		printf("Stereo audio is not supported\n");
		return NULL;
	}
	*/
	// Cancel
	//fread(&channel, sizeof(short int), 1, fp); 
	fread(forIntNumber, sizeof(char), 2, fp);
	channel = forIntNumber[0];
	//printf("\nChannel: %d\n", channel);

	// Sampling frequency
	fread(forIntNumber, sizeof(char), 4, fp);
	*fs = 0;
	for(int i = 3;i >= 0;i--)
	{
	  *fs = (*fs << 8) + forIntNumber[i];
	}
	// Quantization bit rate
	fseek(fp, 6, SEEK_CUR); // Skip 6 bytes
	fread(forIntNumber, sizeof(char), 2, fp);
	*Nbit = forIntNumber[0];
	// Header
	int dummy;
	fread(dataCheck, sizeof(char), 4, fp); // "data"
	while(0 != strcmp(dataCheck,"data"))
	{
		fread(&dummy, sizeof(char), 4, fp);
		fseek(fp, dummy, SEEK_CUR); // Skip irrelevant chunks
		fread(dataCheck, sizeof(char), 4, fp); // "data"
//		fclose(fp);
//		printf("Invalid header data\n");
//		return NULL;
	}
	// Number of sample points
	fread(forIntNumber, sizeof(char), 4, fp); // "data"
	*waveLength = 0;
	for(int i = 3;i >= 0;i--)
	{
	  *waveLength = (*waveLength << 8) + forIntNumber[i];
	}
	*waveLength /= (*Nbit/8 * channel);

	if(*endbr < 0) // If negative, distance from offset
	{
		*endbr = (*waveLength * 1000 / *fs) - (*offset-*endbr);
	}

	int st, ed;
	st = max(0, min(*waveLength-1, (int)((*offset-100) * *fs / 1000))); 
	ed = max(0, min(*waveLength-1, *waveLength - (int)(max(0, *endbr - 100) * *fs / 1000)));
	*endbr = (ed*1000 / *fs) - ((*waveLength * 1000 / *fs) - *endbr);
	*offset = *offset - (st*1000 / *fs);
	*waveLength = (ed - st + 1);

	// Extract waveform
	waveForm = (double *)malloc(sizeof(double) * *waveLength);
	if(waveForm == NULL) return NULL;

	quantizationByte = *Nbit/8;
	zeroLine = pow(2.0,*Nbit-1);
//	for(int i = 0;i < *waveLength;i++)

	fseek(fp, st * quantizationByte * channel, SEEK_CUR);  // Skip to starting position

	unsigned char *wavbuff;
	wavbuff = (unsigned char *) malloc(sizeof(char) * *waveLength * quantizationByte * channel);
	fread(wavbuff, sizeof(char), *waveLength * quantizationByte * channel, fp); // Load into memory
	int seekindex;

	for(int i = 0;i < *waveLength;i++)
	{
		seekindex = i * quantizationByte * channel;
		signBias = 0.0;
		tmp = 0.0;
		// Check sign
		if(wavbuff[seekindex + quantizationByte-1] >= 128)
		{
			signBias = pow(2.0,*Nbit-1);
			wavbuff[seekindex + quantizationByte-1] = wavbuff[seekindex + quantizationByte-1] & 0x7F;
		}
		// Loading data
		for(int j = quantizationByte-1;j >= 0;j--)
		{
			tmp = tmp*256.0 + (double)(wavbuff[seekindex + j]);
		}
		waveForm[i] = (double)((tmp - signBias) / zeroLine);

	}
	// Success!
	free(wavbuff);
	fclose(fp);
	return waveForm;
}

