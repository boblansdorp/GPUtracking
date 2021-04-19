//This will find the new center point using a Fourier Transform Trick
//Cross Correlation
//Using FFTW, doing a 1D Real Array Transform, there are some particular things that happens so read the info here:
//http://www.fftw.org/doc/One_002dDimensional-DFTs-of-Real-Data.html

#include "finddelta.h"

void FindAvgProfileCenter (int cross, int x, int y, double *AvgXProfile, double *AvgYProfile, double *xout, double *yout)
{
	//Initial Variables
	double *inX, *inY;
	complex<double> *outX, *outY;
	fftw_plan pX1, pX2, pY1, pY2;
	int Xmax=0, Ymax=0;

	//This is the input, notice size is Cross
	inX = (double*) fftw_malloc(sizeof(double) * cross);
	inY = (double*) fftw_malloc(sizeof(double) * cross);

	//This is the Complex output of the FFT, notice size is cross/2+1
	outX = (complex<double>*) fftw_malloc(sizeof(fftw_complex) * (cross/2+1));
	outY = (complex<double>*) fftw_malloc(sizeof(fftw_complex) * (cross/2+1));

	//This establishes a plan for the FFTW specifically for 1d real arrays forwards and backwards, for X and Y
	pX1 = fftw_plan_dft_r2c_1d(cross, inX, (fftw_complex*)outX, 0);
    pX2 = fftw_plan_dft_c2r_1d(cross, (fftw_complex*)outX, inX, 0);
	pY1 = fftw_plan_dft_r2c_1d(cross, inY, (fftw_complex*)outY, 0);
    pY2 = fftw_plan_dft_c2r_1d(cross, (fftw_complex*)outY, inY, 0);

	//Assign our current arrays to these
	for(int i=0; i < cross; i++)
	{
		inX[i] = AvgXProfile[i];
		inY[i] = AvgYProfile[i];
	}

	//This executes the FFT
	fftw_execute(pX1);
	fftw_execute(pY1);
	
	//THERE SHOULD BE A COSINE BANDPASS MULTIPLICATION HERE BUT RESULTS INDICATE NOT NEEDED
	//This could be added for higher accuracy?

	//Here we are doing the Cross Correlation Trick
	for(int i =0; i < (cross/2+1); i++)
	{
		outX[i] *= outX[i];
		outY[i] *= outY[i];
	}
	
	fftw_execute(pX2);
	fftw_execute(pY2);

	//Rotate Array so Max point lies in center-ish
	rotate(inX, inX+(cross/2), inX+cross);
	rotate(inY, inY+(cross/2), inY+cross);

	//Normalize the function to match Old VI, all values off by the factor of cross
	//"You should also know that we compute an unnormalized transform" - FFTW Website
	//ALSO ADDED TO THIS LOOP
	//Find the max point (Can this be made smarter? for example we know the max will always be near the center? Do we know this?)
	for(int i = 0; i < cross; i++)
	{
		//we can probably remove this normilization in the end since finding delta is not dependant on this
		AvgXProfile[i] = inX[i] / cross;
		AvgYProfile[i] = inY[i] / cross;

		//Store the index of the max value in Xmax and Ymax
		if(inX[Xmax] < inX[i])
			Xmax = i;
		if(inY[Ymax] < inY[i])
			Ymax = i;

	}

	FindDelta(cross, AvgXProfile, AvgYProfile, Xmax, Ymax, xout, yout);
	
	//RemoveDiscretization(cross, outX, outY);

	*xout = (*xout - (cross/2))/2;
	*yout = (*yout - (cross/2))/2;

	//add delta to original value
	*xout += x;
	*yout += y;

	fftw_destroy_plan(pX1);
	fftw_destroy_plan(pX2);
	fftw_destroy_plan(pY1);
	fftw_destroy_plan(pY2);
	fftw_free(outX); 
	fftw_free(inX);
	fftw_free(outY); 
	fftw_free(inY);
}