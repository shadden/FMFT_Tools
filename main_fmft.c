#include<stdio.h>
#include<math.h>
#include "nrutil.h"

#define PI 3.14159265358979

#define MIN_F     -1.e9   /* minimum frequency (arcsec/yr), 
			    frequencies smaller than MIN_F are ignored */ 
#define MAX_F     1.e9    /* maximum frequency */ 
#define FMFT_FLAG 2      /* the variant of frequency analysis, see fmft.c */   
#define NFREQ     10     /* number of computed frequencies */
#define DATA_SEP  0.08214   /* data spacing (years) */
#define NDATA     8192   /* number of used data, must be a power of 2 and
			    equal or smaller than the number of rows in
			    the input file */ 

int main(void){
  long i;  
  double **output, **input;
  double time, freq, amp, phase;
  double minfreq, maxfreq;

  minfreq = MIN_F /180./3600. * PI *DATA_SEP;
  maxfreq = MAX_F /180./3600. * PI *DATA_SEP;

  input = dmatrix(1,2,1,NDATA);
  output = dmatrix(1,3*FMFT_FLAG,1,NFREQ);

  i=1;
  while(i<=NDATA && scanf("%lg %lg %lg",&time, 
	      &input[1][i],&input[2][i]) != EOF) i++;

  if(i <= NDATA){
    printf("error on input ...\n");
    return 0;
  }

  if(fmft(output, NFREQ, minfreq, maxfreq, FMFT_FLAG, 
	  input, NDATA) == 1) 
 
   for(i=1;i<=NFREQ;i++){
      freq = output[3*FMFT_FLAG-2][i] *180.*3600./PI / DATA_SEP; 
                                             /* arsec per unit of DATA_SEP */
      amp = output[3*FMFT_FLAG-1][i];
      phase = output[3*FMFT_FLAG][i] *180./PI;  /* degrees */
      if(phase < 0) phase += 360.;
	
      printf("%.15lg %.15lg %.15lg\n", freq, amp, phase);

    }
  
  return 1;
}









