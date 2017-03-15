#include <stdio.h>
int fmft(double **output, int nfreq, double minfreq, double maxfreq, int flag, 
	 double **input, long ndata);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);

int fmftWrap(double my_out[][3], int nfreq, double minfreq, double maxfreq, int flag, 
	double my_in[][3], long ndata)
{
	double **output; double **input;
	double dt;
	int i,err;
	
	
	input = dmatrix(1,2,1,ndata);
	output = dmatrix(1,3*flag,1,nfreq);
	
	i =0;
	while(i<ndata){
		input[1][i] = my_in[i][1];
		input[2][i] = my_in[i][2];
		i++;
	}
	
	err = fmft(output, nfreq, minfreq, maxfreq, flag, input, ndata);
	
	
	dt = my_in[1][0] - my_in[0][0];
	

	for(i=0;i<nfreq; i++){
		
		// Frequency
		my_out[i][0] = output[3*flag-2][i+1] / dt;
		// Amplitude
		my_out[i][1] = output[3*flag-1][i+1];
		// Phase 
		my_out[i][2] = output[3*flag][i+1];
	}

	
	free_dmatrix(input,1,2,1,ndata);
	free_dmatrix(output,1,3*flag,1,nfreq);
	
	return err;
}

double array2D(double array[][2], int num_entries)
{
    int i;
    double sum = 0.0;
    for (i = 0; i < num_entries; i++)
    {
        sum += array[i][0] + array[i][1];
    }
    return sum ;
}