#define SIZE 5 //this is the number of nearest neighbouring slices in the calibration stack that are used for fitting... 5 works

//this finds the least squares fit of the prepped I(r) to find the index of the closest slice
__global__ void FitPrepandPhase (int cross, int *forgetradius, int numofimages, int numofbeads, int *RealRow, int *RealCol, int *AmpRow, int *AmpCol, int *CorkRow, int *CorkCol, cufftComplex* radinprofile, float *Real, float *amplitudes, cufftComplex *corkscrews, float *zout, bool* beadstatus)
{
	int ImageID = threadIdx.x + blockIdx.x * blockDim.x;
	int BeadID = threadIdx.y + blockIdx.y * blockDim.y;
if ((ImageID < numofimages) && (BeadID < numofbeads) && beadstatus[BeadID])
{
	//We have to run functions like this to take into account the variable array sizes
	int TotalRealDisp=0, TotalCorkDisp=0, TotalAmpDisp=0;
	for (int i=0; i < BeadID; i++)
	{
		TotalRealDisp += RealRow[i]*RealCol[i];
		TotalCorkDisp += CorkRow[i] * CorkCol[i];
		TotalAmpDisp += AmpRow[i] * AmpCol[i];
	}

	//We create Prep and Corkscrew Arrays
	float *prep = new float[((cross/2)-forgetradius[BeadID])];
	//float prep[((cross/2)-forgetradius[BeadID])];
	cufftComplex *corkscrew = new cufftComplex[((cross/2)-forgetradius[BeadID])];
	for(int i=0; i < ((cross/2)-forgetradius[BeadID]); i++)
	{
		prep[i] = radinprofile[(ImageID*numofbeads*(cross-1))+(BeadID*(cross-1))+ (((cross/2)-forgetradius[BeadID])-1)-i].x/(cross-1);
		corkscrew[i].x = radinprofile[(ImageID*numofbeads*(cross-1))+(BeadID*(cross-1)) + (((cross/2)-forgetradius[BeadID])-1)-i].x/(cross-1);
		corkscrew[i].y = -(radinprofile[(ImageID*numofbeads*(cross-1))+(BeadID*(cross-1)) + (((cross/2)-forgetradius[BeadID])-1)-i].y/(cross-1));
	}
	
	//We use least square to find the closest fit of prep on the premade Real Matrix
	int index=0;
	float *Arr = new float[RealRow[BeadID]];
	//float Arr[RealRow[BeadID]];
	for (int i=0; i < RealRow[BeadID]; i++)
	{
		Arr[i]=0;
		for (int j=0; j < RealCol[BeadID]; j++)
		{
			//Arr[i] += pow((prep[j] - Real[(TotalRealDisp) + ((i * RealCol[BeadID]) + j)]), 2);
			Arr[i] += (prep[j] - Real[(TotalRealDisp) + ((i * RealCol[BeadID]) + j)])*(prep[j] - Real[(TotalRealDisp) + ((i * RealCol[BeadID]) + j)]);
		}
		if(Arr[i] < Arr[index])
			index = i;
	}
	delete[] Arr;
	delete[] prep;

	//Using Complex numbers:
	//complex in polar form is x = r * e(i* theta)
	
	//Use the phase information from the closes fit above
	//int size = 5;
	float phases[SIZE]; //float *phases = new float[size];
	int indexstart = (int)(index - (SIZE/2));
	if (indexstart<0)
		indexstart=0; //prevent negative values of array index later (memory read error)

	int length = ((cross/2)-forgetradius[BeadID]);
	float temp;
	for(int i=0; i < SIZE; i++)
	{
		float sum1=0, sum2=0;
		for(int j=0; j < length; j++)
		{
			float a=0, b=0, c=0, d=0, e=0;
			//This is written this way to be easier to read, no reason to make all these variables.
			//We are manipulating complex numbers in polar form.
			a = corkscrew[j].x;
			b = corkscrew[j].y;
			c = corkscrews[(TotalCorkDisp) + (((i+indexstart) * CorkCol[BeadID]) + j)].x;
			d = corkscrews[(TotalCorkDisp) + (((i+indexstart) * CorkCol[BeadID]) + j)].y;
			e = amplitudes[(TotalAmpDisp) + (((i+indexstart) * AmpCol[BeadID]) + j)];
			temp = sqrt(a*a + b*b) * e;
			sum1 += temp;
			//Changed this from atan to atan2 because atan2 has a larger range and will give a more accurate output.
			sum2 += temp * atan2((b*c - a*d),(a*c + b*d));
		}
		phases[i]=sum2/sum1;
	}

	delete[] corkscrew;

	//This takes a descrete value for the best fit cal image index, and then uses a quadratic fit
	//on the phase information to find a non integer more accurate value for the best fit z height

	float data[SIZE*2]; //float *data = new float[size*2];
	float solution[3] = {0};
	for(int i=0; i < SIZE; i++)
	{
		//note that phases is the domain
		data[i*2+0] = phases[i];
		data[i*2+1] = (i - floor((float)SIZE/2));
	}

	QuadFit(SIZE, data, solution);
	zout[(ImageID*numofbeads)+BeadID] = (index + solution[2]);
	//delete[] phases;
	//delete[] data;
}
else if((ImageID < numofimages) && (BeadID < numofbeads) && !beadstatus[BeadID])
{
	zout[(ImageID*numofbeads)+BeadID] = 0;
}
}

//occasionally, 