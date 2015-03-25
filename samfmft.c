#include<stdio.h>
#include<math.h>
#include "nrutil.h"
#include <assert.h>

#define PI 3.14159265358979

#define MIN_F     -20.e0   /* minimum frequency (arcsec/yr), 
			    frequencies smaller than MIN_F are ignored */ 
#define MAX_F     20.e0    /* maximum frequency */ 
#define FMFT_FLAG 2      /* the variant of frequency analysis, see fmft.c */   
#define NFREQ     4     /* number of computed frequencies */
#define DATA_SEP  1.   /* data spacing (years) */
#define NDATA     2048  /* number of used data, must be a power of 2 and
			    equal or smaller than the number of rows in
			    the input file */ 

int main(int argc, char *argv[]){
  long i;  
  double **output, **input;
  double time, freq, amp, phase;
  double minfreq, maxfreq;
  double dsep;

  assert(argc == 2 );
  sscanf(argv[1],"%lg",&dsep);
  printf("dsep set to: %lg \n",dsep);
  fflush(stdout);

  minfreq = MIN_F * dsep;
  maxfreq = MAX_F * dsep;

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
      freq = output[3*FMFT_FLAG-2][i] / dsep  /* radians per  */ ; 
                                             
      amp = output[3*FMFT_FLAG-1][i];
      phase = output[3*FMFT_FLAG][i] * 180./PI;  /* degrees */
      if(phase < 0) phase += 360.;
	
      printf("%.15lg %.15lg %.15lg\n", freq, amp, phase);

    }
  
  return 1;
}
