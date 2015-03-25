#include "mathlink.h"

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define PI 3.14159265358979

#define MIN_F     -20.e0   /* minimum frequency (arcsec/yr), 
			    frequencies smaller than MIN_F are ignored */ 
#define MAX_F     20.e0    /* maximum frequency */ 
#define FMFT_FLAG 2      /* the variant of frequency analysis, see fmft.c */   
#define DATA_SEP  1.   /* data spacing (years) */


int fmft(double **output, int nfreq, double minfreq, double maxfreq, int flag, 
	 double **input, long ndata);

double **dmatrix(long nrl, long nrh, long ncl, long nch);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);

void throw_error(char err_msg[200]){
		MLEvaluate(stdlink, err_msg); 
 		MLNextPacket(stdlink);
		MLNewPacket(stdlink);
		MLPutSymbol(stdlink,"$Failed");
		return;
	}
	
void ML_Attempt(int func){
	if (! func ){
		char err_msg[200];
		MLClearError(stdlink);
		sprintf(err_msg, "Message[%s,\"%.76s\"]","FrequencyModifiedFourierTransform::mlink", MLErrorMessage(stdlink));
		throw_error(err_msg);
		return;
	}
	return;
}

void FrequencyModifiedFourierTransform(int nfreq){
	// Variables for Mathematica interface
	double *data;
	int *dimensions;
	int depth;
	char **heads;
	
	
	// Variables for FMFT
	long i, length_pow2;
	int length;
	double **output, **input;
	double time, freq, amp, phase;
	double minfreq, maxfreq;
	double dsep;
	
	
	
	// Read in mathematica data array
	if( ! MLGetReal64Array(stdlink,&data,&dimensions,&heads,&depth) ){
		char err_msg[300];
		sprintf(err_msg, "%s\"%.76s\"%s","Message[FrequencyModifiedFourierTransform::mlink,",MLErrorMessage(stdlink),"]");
	 	throw_error(err_msg);
		return;
	 }
	 // ... ensure that input array has the appropriate dimensions
	 if( depth != 2 ){
		char err_msg[200];
		sprintf(err_msg, "%s","Message[FrequencyModifiedFourierTransform::BadDimensions]");
	 	throw_error(err_msg);	 
		return;
	}
	
	// Get closest power of 2 to length of input data
	length = dimensions[0];
	length_pow2 = (long) pow(2. , floor( log2( (float) length ) ) );
	
	// Compute the time separation between successive data points and
	//  calculate the minimum and maximum frequency based in this separation
	dsep = data[0+1*dimensions[1]] - data[0+0*dimensions[1]];
	// Set min/max frequency a little slower than nyquist frequency, which should be 0.5 (?)
	minfreq = -0.3; // MIN_F * dsep;
  	maxfreq = 0.3;  // MAX_F * dsep;
	
	input = dmatrix( 1,2,1,length_pow2 );
	output = dmatrix(1,3*FMFT_FLAG,1,nfreq);
	
	i=0 ;
	while(i< length_pow2){
		input[1][i] = data[1+i*dimensions[1]];
		input[2][i] = data[2+i*dimensions[1]];
		i++;
	}
	
	if(fmft(output, nfreq, minfreq, maxfreq, FMFT_FLAG, input, length_pow2) != 1){
		char err_msg[200];
		sprintf(err_msg, "%s","Message[FrequencyModifiedFourierTransform::FMFTFailure]");
	 	throw_error(err_msg);
		return;
	}
	// Convert frequency values to proper units based in data spacing

	
	ML_Attempt( MLPutFunction(stdlink,"List",nfreq) );
	for(i=1;i<=nfreq; i++){
		double outarr[3];
		// Frequency
		outarr[0] = output[3*FMFT_FLAG-2][i] / dsep;
		// Amplitude
		outarr[1] = output[3*FMFT_FLAG-1][i];
		// Phase 
		outarr[2] = output[3*FMFT_FLAG][i];
		// Send output info to Mathematica
		ML_Attempt( MLPutReal64List(stdlink,outarr,3) );
	}
	
	// Free allocated variables
	free_dmatrix(input,1,2,1,length_pow2);
	free_dmatrix(output,1,3*FMFT_FLAG,1,nfreq);
	// Release memory for Mathematica input array
	MLReleaseReal64Array(stdlink,data,dimensions,heads,depth);
	
	return ;
}

int main(int argc, char *argv[]){
	return MLMain(argc,argv);
}
