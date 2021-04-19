//This will find delta using quadfit and maybe some fourier transformas and manipulations in the future

#include "quadfit.h"

void FindDelta (int cross, double *AvgXProfile, double *AvgYProfile, int Xmax, int Ymax, double *xout, double *yout)
{
	//solution should come from quadfit, and bein the form ax^2 + bx + c = y
	//a = solution [0], b = solution[1], c = solution[2]
	//to find a max, we want 2ax + b = 0 => x = -b/2a which is our new starting location
	

	double solutionX[3] = {0};
	double solutionY[3] = {0};

	//determines size of subarray
	int size = 5;
	double **subinX, **subinY;

	//These will be our subarrays to pass into quadfit
	subinX = new double*[size];
	subinY = new double*[size];
	for(int i=0; i < size; i++)
	{
		subinX[i] = new double[2];
		subinY[i] = new double[2];

		subinX[i][0] = (int)((Xmax - (size/2))+i);
		subinX[i][1] = AvgXProfile[(int)((Xmax - (size/2))+i)];
		subinY[i][0] = (int)((Ymax - (size/2))+i);
		subinY[i][1] = AvgYProfile[(int)((Ymax - (size/2))+i)];
	}

	QuadFit(size, subinX, solutionX);
	QuadFit(size, subinY, solutionY);


	*xout = -(solutionX[1]/(2*solutionX[0]));
	*yout = -(solutionY[1]/(2*solutionY[0]));

	for(int i=0; i < size; i++)
	{
		delete[] subinX[i];
		delete[] subinY[i];
	}

	delete[] subinX;
	delete[] subinY;
}