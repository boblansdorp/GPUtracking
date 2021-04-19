//This function removes descretization problems that arise in the XY Tracking of the beads when trying to resolve the position to 1/100th of a pixel.
//This function is an original idea by Vincent Croquette from the ABCD Biophysics Lab
//It was originally impimented in LabView and modified for CUDA for this particular program.

//Vincent Croquette
//vincent@lps.ens.fr
//http://pimprenelle.lps.ens.fr/biolps/user/vincent-croquette

__global__ void OriginalFit(int numofimages, int numofbeads, float *delta, float *fit)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
while (tid < numofimages*numofbeads)
{
	//This is one method for calulating the nearest integer value of some float.
	//However we should probaly repalce this with a nonbiased version.
	//In labview, the implimentation they use is round to nearest even integer.
	fit[tid] = floor(delta[tid] + .5);
	tid += blockDim.x * gridDim.x;
}
}

__global__ void FitFix(int numofimages, int numofbeads, float *out, float *fit)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
while (tid < numofimages*numofbeads)
{
	out[tid] -= fit[tid];
	tid += blockDim.x * gridDim.x;
}
}

__global__ void PhaseTrick(int cross, int numofimages, int numofbeads, float *delta, float *fit, cufftComplex *AvgProfileOut)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;

	float rad, ang, a, fitdec;
	//Here we use the displacement found from FindDelta, and the orignal fit to shift the phase of the AvgProfile.
	while (tid < cross * numofimages * numofbeads)
	{
		//This is the difference between the nearest integer value fit and the displacement, giving us just the decimal value.
		fitdec = delta[tid/cross] - fit[tid/cross];
		if ((tid%cross) > (cross/2))
		{
			a = ((tid%cross) - cross);
		} else
		{
			a = (tid%cross);
		}
		//Conversion from Complex to Polar
		rad = sqrt(AvgProfileOut[tid].x * AvgProfileOut[tid].x + AvgProfileOut[tid].y * AvgProfileOut[tid].y);
		ang = atan2(AvgProfileOut[tid].y, AvgProfileOut[tid].x);
		//Shift the phase
		ang += ((2*PI/cross)*fitdec*a);
		//Convert back to Complex
		AvgProfileOut[tid].x = rad * cos(ang);
		AvgProfileOut[tid].y = rad * sin(ang);
		
		tid += blockDim.x * gridDim.x;
	}
}

void RemoveDiscretization(int maxthreadsperblock, dim3 numBlocks, dim3 threadsPerBlock, int cross, int NUMOFIMAGES, int NUMOFBEADS, cufftHandle plan, cufftComplex *devAvgProfileIn, cufftComplex *devAvgProfileOut, float* devfit, float*devdelta, float* devout, bool* devbeadstatus)
{
	//Using the displacement found in FindDelta, we shift the phase of the AvgProfile and then recaclulate the displalcement and add that as part of our XY position.
	PhaseTrick<<<(1+(cross*NUMOFIMAGES*NUMOFBEADS)/maxthreadsperblock), maxthreadsperblock>>>(cross, NUMOFIMAGES, NUMOFBEADS, devdelta, devfit, devAvgProfileOut);

	cufftExecC2C(plan, devAvgProfileOut, devAvgProfileIn, CUFFT_INVERSE);

	Rotate<<<numBlocks,threadsPerBlock>>>(NUMOFIMAGES, NUMOFBEADS, cross, devAvgProfileIn);

	FindDelta<<<numBlocks,threadsPerBlock>>>(cross, NUMOFIMAGES, NUMOFBEADS, devAvgProfileIn, devdelta, devout, devbeadstatus);
	
	//We need to make sure to remove 'the original fit' which is the nearest integer value to the original displacement to properly convey the location of the bead.
	FitFix<<<(1+(NUMOFIMAGES*NUMOFBEADS)/maxthreadsperblock), maxthreadsperblock>>>(NUMOFIMAGES, NUMOFBEADS, devout, devfit);
}