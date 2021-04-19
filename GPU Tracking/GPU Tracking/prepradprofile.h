//cosine bandpass application
__global__ void CosBandPass (int cross, int numofimages, int numofbeads, cufftComplex *radinprofile, float *cosband, int *CosBandLength)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	while (tid < numofimages*numofbeads*(cross-1))
	{
		int TotalCosBandLength = 0;
		int beadnum = floor((float)(tid%(numofbeads*(cross-1)))/(cross-1));

		for(int i=0; i < beadnum; i++)
			TotalCosBandLength += CosBandLength[i];

		radinprofile[tid].x = radinprofile[tid].x * cosband[TotalCosBandLength+tid%(cross-1)];
		radinprofile[tid].y = radinprofile[tid].y * cosband[TotalCosBandLength+tid%(cross-1)];	
		tid += blockDim.x * gridDim.x;
	}
}

//This program preps the radial profile so that calculations can be made on it
__host__ void PrepRadInProf (int maxthreadsperblock, int cross, int numofimages, int numofbeads, cufftComplex *radinprofile, float *cosband, int *CosBandLength)
{
	cufftHandle plan;
	cufftPlan1d (&plan, (cross-1), CUFFT_C2C, numofimages*numofbeads);

	//forward transform into hilbert space
	cufftExecC2C(plan, radinprofile, radinprofile, CUFFT_FORWARD);

	//Apply the cosine bandpass to the entire array
	CosBandPass<<<(1+((cross-1)*numofimages*numofbeads)/maxthreadsperblock), maxthreadsperblock>>>(cross, numofimages, numofbeads, radinprofile, cosband, CosBandLength);
	
	//reverse transform back
	cufftExecC2C(plan, radinprofile, radinprofile, CUFFT_INVERSE);

	cufftDestroy ( plan ) ;

}