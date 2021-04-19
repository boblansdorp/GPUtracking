DllExport int GPUTracking (int threads, int DIM1, int DIM2, int DIM3, int DIM4, Array1dIntHandle x, Array1dIntHandle y, int cross, int crossthickness, ImageHandle image, Array2dHandle xout, Array2dHandle yout, Array2dHandle zout, ArrayClusterHandle calibration, BoolArrayHandle BeadStatus, char* text, Array1dHandle test)
{
	cudaError_t error = cudaSuccess;

	int NUMOFBEADS = (*x)->length;
	int NUMOFIMAGES = (*image)->length[0];
	int NUMROWS = (*image)->length[1];
	int NUMCOLS = (*image)->length[2];

	//TODO: Create a function that automatically optimizes the parmeters below based off the CUDA Deivice Properties
	//Must take into account the Number of Beads, Threads, Max Number of Threads, Max Registers per Block, etc...
	int maxthreadsperblock = threads;
	if (threads == 0)
	{
		cudaDeviceProp deviceProp;
		cudaGetDeviceProperties(&deviceProp, 0);
		maxthreadsperblock = deviceProp.maxThreadsPerBlock;
	//	maxthreadsperblock/=4;
	} else
	{
		maxthreadsperblock = threads;
	}
	(*test)->val[0] = (float)maxthreadsperblock;

	//Depending on your GPU, you might have a limited number of threads per block you can allocate.
	//Take a look at the function above that outputs the Max Threads Per Block.
	//You will also have to adjust this function for the experiment you are performing.
	//The first dimension is directly related to the number of images, and the second dimension is related to number of beads.
	//If you are performing an experiment with only 2 beads, you should not have the second dimension larger than 2 or you will waste threads.
	//Similarly with number of images. Remember you must have: DIM1 * DIM2 < Max Number of Threads
	dim3 threadsPerBlock( DIM1, DIM2  );
	//Make enough blocks such that every single bead has its own thread
	dim3 numBlocks( (1 + (NUMOFIMAGES-1)/threadsPerBlock.x), (1 + (NUMOFBEADS-1)/threadsPerBlock.y));

	//Calculate Radial Profile is independently optimized to perform better, since it is the slowest of all the functions.
	//Notice this now uses 3 dimensions, and the 3rd dimension should ALWAYS remain "cross"
	//Similarly with the above optimization, DIM1 and DIM2 are directly related to Num of Images and Threads.
	//However you must now have DIM1 * DIM2 * Cross < MaxNumberOfThreads
	//threadsRadial(cross, beads, images)
	//dim3 threadsRadial(cross, DIM3, DIM4);


	dim3 threadsRadial(cross, DIM3, DIM4);


	dim3 blocksRadial( 1, (1 + (NUMOFIMAGES-1)/threadsRadial.y), (1 + (NUMOFBEADS-1)/threadsRadial.z));
	
	
	dim3 threadsfitprep(1, NUMOFBEADS,1);
	dim3 blocksfitprep( (1 + (NUMOFIMAGES-1)/threadsPerBlock.x), (1 + (NUMOFBEADS-1)/threadsPerBlock.x),1);

	
	

//Start XY Tracking
//All XY Tracking could be parellized to do the X and Y tracking process at the same time since they are compeltely independent.
	unsigned char *devimage;
	float *devxout, *devyout, *devxdelta, *devydelta, *devxfit, *devyfit;
	bool *devbeadstatus;
	cufftComplex *devAvgXProfileIn, *devAvgYProfileIn, *devAvgXProfileOut, *devAvgYProfileOut;
	int *devxin, *devyin;

	//Allocate Memory	
	cudaMalloc( (void**)&devimage, sizeof(unsigned char) * NUMOFIMAGES*NUMCOLS*NUMROWS);
	cudaMalloc( (void**)&devbeadstatus, sizeof(bool) * NUMOFBEADS);
	cudaMalloc( (void**)&devxin, sizeof(int) * NUMOFBEADS);
	cudaMalloc( (void**)&devyin, sizeof(int) * NUMOFBEADS);

	cudaMalloc( (void**)&devxout, sizeof(float)* NUMOFIMAGES* NUMOFBEADS );
	cudaMalloc( (void**)&devyout, sizeof(float)* NUMOFIMAGES* NUMOFBEADS );

	cudaMalloc( (void**)&devAvgXProfileIn, sizeof(cufftComplex) * cross * NUMOFIMAGES * NUMOFBEADS);
	cudaMalloc( (void**)&devAvgYProfileIn, sizeof(cufftComplex) * cross * NUMOFIMAGES * NUMOFBEADS);
	cudaMalloc( (void**)&devAvgXProfileOut, sizeof(cufftComplex) * cross * NUMOFIMAGES * NUMOFBEADS);
	cudaMalloc( (void**)&devAvgYProfileOut, sizeof(cufftComplex) * cross * NUMOFIMAGES * NUMOFBEADS);

	cudaMalloc( (void**)&devxfit, sizeof(float)* NUMOFIMAGES* NUMOFBEADS );
	cudaMalloc( (void**)&devyfit, sizeof(float)* NUMOFIMAGES* NUMOFBEADS );

	cudaMalloc( (void**)&devxdelta, sizeof(float)* NUMOFIMAGES* NUMOFBEADS );
	cudaMalloc( (void**)&devydelta, sizeof(float)* NUMOFIMAGES* NUMOFBEADS );

	//Copy Values from LabView
	cudaMemcpy( devimage, (*image)->val, sizeof(unsigned char) * NUMOFIMAGES*NUMCOLS*NUMROWS, cudaMemcpyHostToDevice);
	cudaMemcpy( devbeadstatus, (*BeadStatus)->val, sizeof(bool) * NUMOFBEADS, cudaMemcpyHostToDevice);
	cudaMemcpy( devxin, (*x)->val, NUMOFBEADS * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy( devyin, (*y)->val, NUMOFBEADS * sizeof(int), cudaMemcpyHostToDevice);

	//Create and Prep the X and Y Cross Profile
	AvgXCrossAndPrep<<<numBlocks,threadsPerBlock>>>(devxin, devyin, cross, crossthickness, NUMOFIMAGES, NUMOFBEADS, NUMROWS, NUMCOLS, devimage, devAvgXProfileIn, devbeadstatus, devxout);
	AvgYCrossAndPrep<<<numBlocks,threadsPerBlock>>>(devxin, devyin, cross, crossthickness, NUMOFIMAGES, NUMOFBEADS, NUMROWS, NUMCOLS, devimage, devAvgYProfileIn, devbeadstatus, devyout);
	
	//Create a plan for the FFT. Do NUMOFIMAGES*NUMOFBEADS number of FFTs each of length cross.
	cufftHandle plan;
	cufftPlan1d(&plan, cross, CUFFT_C2C, NUMOFIMAGES*NUMOFBEADS);
	
	//Execute the FFT
	cufftExecC2C(plan, devAvgXProfileIn, devAvgXProfileOut, CUFFT_FORWARD);
	cufftExecC2C(plan, devAvgYProfileIn, devAvgYProfileOut, CUFFT_FORWARD);

	//Note that this function uses as many threads as possible to perform this task, not the grid defined above.
	//Square all values after the FFT, amd KEEP THESE VALUES
	SquareVector<<<(1+(cross*NUMOFIMAGES*NUMOFBEADS)/maxthreadsperblock), maxthreadsperblock>>>(cross, NUMOFIMAGES, NUMOFBEADS, devAvgXProfileOut);
	SquareVector<<<(1+(cross*NUMOFIMAGES*NUMOFBEADS)/maxthreadsperblock), maxthreadsperblock>>>(cross, NUMOFIMAGES, NUMOFBEADS, devAvgYProfileOut);
	
	//Reverse transform back into AvgProfileIn, while maintaining the squared output vector.
	cufftExecC2C(plan, devAvgXProfileOut, devAvgXProfileIn, CUFFT_INVERSE);
	cufftExecC2C(plan, devAvgYProfileOut, devAvgYProfileIn, CUFFT_INVERSE);

	//Rotate all the profiles so that the PEAK is at the center rather than the start of the array
	Rotate<<<numBlocks,threadsPerBlock>>>(NUMOFIMAGES, NUMOFBEADS, cross, devAvgXProfileIn);
	Rotate<<<numBlocks,threadsPerBlock>>>(NUMOFIMAGES, NUMOFBEADS, cross, devAvgYProfileIn);

	//Find the displacement of the center of the bead to our original guess location. Also save the displacement to be used for further corrections.
	FindDelta<<<numBlocks,threadsPerBlock>>>(cross, NUMOFIMAGES, NUMOFBEADS, devAvgXProfileIn, devxdelta, devxout, devbeadstatus);
	FindDelta<<<numBlocks,threadsPerBlock>>>(cross, NUMOFIMAGES, NUMOFBEADS, devAvgYProfileIn, devydelta, devyout, devbeadstatus);

	//Take the original displacement (delta) that we find, and round it to the nearest integer to fix descretization.
	OriginalFit<<<(1+(NUMOFIMAGES*NUMOFBEADS)/maxthreadsperblock), maxthreadsperblock>>>(NUMOFIMAGES, NUMOFBEADS, devxdelta, devxfit);
	OriginalFit<<<(1+(NUMOFIMAGES*NUMOFBEADS)/maxthreadsperblock), maxthreadsperblock>>>(NUMOFIMAGES, NUMOFBEADS, devydelta, devyfit);

	//This loop/function helps fix Discretization issues that arise when trying to resolve the bead position to 1/100th of a pixel. We run it twice.
	
	for (int i=0; i < 2; i++)
	{
		RemoveDiscretization(maxthreadsperblock, numBlocks, threadsPerBlock, cross, NUMOFIMAGES, NUMOFBEADS, plan, devAvgXProfileIn, devAvgXProfileOut, devxfit, devxdelta, devxout, devbeadstatus);
		RemoveDiscretization(maxthreadsperblock, numBlocks, threadsPerBlock, cross, NUMOFIMAGES, NUMOFBEADS, plan, devAvgYProfileIn, devAvgYProfileOut, devyfit, devydelta, devyout, devbeadstatus);
	}	

	//We correct our XY position to reflect its position on the Image.
	Out<<<numBlocks,threadsPerBlock>>>(NUMOFIMAGES, NUMOFBEADS, cross, devxin, devxout, devbeadstatus);
	Out<<<numBlocks,threadsPerBlock>>>(NUMOFIMAGES, NUMOFBEADS, cross, devyin, devyout, devbeadstatus);
	
	//Copy the actual XY position back to LabView
	cudaMemcpy( (*xout)->val, devxout, NUMOFIMAGES*NUMOFBEADS*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy( (*yout)->val, devyout, NUMOFIMAGES*NUMOFBEADS*sizeof(float), cudaMemcpyDeviceToHost);

	cufftDestroy(plan);

//End XY Tracking


//Start Z Tracking
	//The following process allows us to use a different Forget Radius for each bead. The size of forget radius affects the length of these other arrays.
	int MaxAmpRow = 0, MaxAmpCol = 0, MaxCorkRow = 0, MaxCorkCol = 0, MaxRealRow = 0, MaxRealCol = 0, MaxCosBandLength = 0;
	for(int i=0; i < NUMOFBEADS; i++)
	{
		if(MaxAmpRow < (*(*calibration)->cluster[i].amplitude)->length[0])
			MaxAmpRow = (*(*calibration)->cluster[i].amplitude)->length[0];
		if(MaxAmpCol < (*(*calibration)->cluster[i].amplitude)->length[1])
			MaxAmpCol = (*(*calibration)->cluster[i].amplitude)->length[1];
		if(MaxCorkRow < (*(*calibration)->cluster[i].corkscrews)->length[0])
			MaxCorkRow = (*(*calibration)->cluster[i].corkscrews)->length[0];
		if(MaxCorkCol < (*(*calibration)->cluster[i].corkscrews)->length[1])
			MaxCorkCol = (*(*calibration)->cluster[i].corkscrews)->length[1];
		if(MaxRealRow < (*(*calibration)->cluster[i].Real)->length[0])
			MaxRealRow = (*(*calibration)->cluster[i].Real)->length[0];
		if(MaxRealCol < (*(*calibration)->cluster[i].Real)->length[1])
			MaxRealCol = (*(*calibration)->cluster[i].Real)->length[1];
		if(MaxCosBandLength < (*(*calibration)->cluster[i].cosband)->length)
			MaxCosBandLength = (*(*calibration)->cluster[i].cosband)->length;
	}

	float *devzout, *devcosband, *devReal, *devamplitudes;
	int *devforgetradius, *devAmpRow, *devAmpCol, *devCorkRow, *devCorkCol, *devRealRow, *devRealCol, *devCosBandLength;
	cufftComplex *devcorkscrews, *devradprofileout;

	//Allocate Memory
	//Notice we OVERALLOCATE memory to take into account the varying size of the arrays.
	cudaMalloc( (void**)&devcosband, sizeof(float) * NUMOFBEADS * MaxCosBandLength);
	cudaMalloc( (void**)&devradprofileout, sizeof(cufftComplex) * NUMOFBEADS * (cross-1) * NUMOFIMAGES);
	cudaMalloc( (void**)&devzout, sizeof(float) * NUMOFBEADS * NUMOFIMAGES);
	cudaMalloc( (void**)&devforgetradius, sizeof(int) * NUMOFBEADS);
	cudaMalloc( (void**)&devAmpRow, sizeof(int) * NUMOFBEADS);
	cudaMalloc( (void**)&devAmpCol, sizeof(int) * NUMOFBEADS);
	cudaMalloc( (void**)&devCorkRow, sizeof(int) * NUMOFBEADS);
	cudaMalloc( (void**)&devCorkCol, sizeof(int) * NUMOFBEADS);
	cudaMalloc( (void**)&devRealRow, sizeof(int) * NUMOFBEADS);
	cudaMalloc( (void**)&devRealCol, sizeof(int) * NUMOFBEADS);
	cudaMalloc( (void**)&devCosBandLength, sizeof(int) * NUMOFBEADS);
	cudaMalloc( (void**)&devcorkscrews, sizeof(cufftComplex) * NUMOFBEADS * MaxCorkRow * MaxCorkCol);
	cudaMalloc( (void**)&devReal, sizeof(float) * NUMOFBEADS * MaxRealRow * MaxRealCol);
	cudaMalloc( (void**)&devamplitudes, sizeof(float) * NUMOFBEADS * MaxAmpRow * MaxAmpCol);

	//The following loop copies values from LabView to the GPU Memory while taking into account the variable lengths of the arrays.
	int CorkDisp=0, RealDisp=0, AmpDisp=0, CosDisp=0;
	for (int i=0; i< NUMOFBEADS; i++)
	{
		int AmpRow = (*(*calibration)->cluster[i].amplitude)->length[0];
		int AmpCol = (*(*calibration)->cluster[i].amplitude)->length[1];
		int CorkRow = (*(*calibration)->cluster[i].corkscrews)->length[0];
		int CorkCol = (*(*calibration)->cluster[i].corkscrews)->length[1];
		int RealRow = (*(*calibration)->cluster[i].Real)->length[0];
		int RealCol = (*(*calibration)->cluster[i].Real)->length[1];
		int CosBandLength = (*(*calibration)->cluster[i].cosband)->length;

		cudaMemcpy( &devforgetradius[i], &(*calibration)->cluster[i].forgetradius, sizeof(int), cudaMemcpyHostToDevice );
		cudaMemcpy( &devAmpRow[i], &(*(*calibration)->cluster[i].amplitude)->length[0], sizeof(int), cudaMemcpyHostToDevice );
		cudaMemcpy( &devAmpCol[i], &(*(*calibration)->cluster[i].amplitude)->length[1], sizeof(int), cudaMemcpyHostToDevice );
		cudaMemcpy( &devCorkRow[i], &(*(*calibration)->cluster[i].corkscrews)->length[0], sizeof(int), cudaMemcpyHostToDevice );
		cudaMemcpy( &devCorkCol[i], &(*(*calibration)->cluster[i].corkscrews)->length[1], sizeof(int), cudaMemcpyHostToDevice );
		cudaMemcpy( &devRealRow[i], &(*(*calibration)->cluster[i].Real)->length[0], sizeof(int), cudaMemcpyHostToDevice );
		cudaMemcpy( &devRealCol[i], &(*(*calibration)->cluster[i].Real)->length[1], sizeof(int), cudaMemcpyHostToDevice );
		cudaMemcpy( &devCosBandLength[i], &(*(*calibration)->cluster[i].cosband)->length, sizeof(int), cudaMemcpyHostToDevice );
		
		cudaMemcpy( &devcorkscrews[CorkDisp], (*(*calibration)->cluster[i].corkscrews)->val, sizeof(cufftComplex)*CorkRow*CorkCol, cudaMemcpyHostToDevice );
		CorkDisp += CorkRow*CorkCol;
		cudaMemcpy( &devReal[RealDisp], (*(*calibration)->cluster[i].Real)->val, sizeof(float) * RealRow * RealCol, cudaMemcpyHostToDevice );
		RealDisp += RealRow * RealCol;
		cudaMemcpy( &devamplitudes[AmpDisp], (*(*calibration)->cluster[i].amplitude)->val, sizeof(float) * AmpRow * AmpCol, cudaMemcpyHostToDevice );
		AmpDisp += AmpRow * AmpCol;
		cudaMemcpy( &devcosband[CosDisp], (*(*calibration)->cluster[i].cosband)->val, sizeof(float) * CosBandLength, cudaMemcpyHostToDevice);
		CosDisp += CosBandLength;
	}

	//Calculate Radial Profile is one of the most important, and most time consuming functions
	//We have attempted to optimize it significantly, and have succeeded in making it run significantly faster, yet it may be possible for futher work to still be done.
	float *devradinprofile;
	cudaMalloc( (void**)&devradinprofile, sizeof(float)* NUMOFIMAGES* NUMOFBEADS*cross);

	//Note that Calculate Radial Profile 1 uses a unique grid to make blocks and threads, this is the function we expect to take alot of time.
	CalculateRadialProfile1<<<blocksRadial,threadsRadial>>>(cross, NUMOFIMAGES, NUMOFBEADS, NUMCOLS, NUMROWS, devxout, devyout, devradprofileout, devimage, devbeadstatus, devradinprofile);
	CalculateRadialProfile2<<<numBlocks,threadsPerBlock>>>(cross, NUMOFIMAGES, NUMOFBEADS, NUMCOLS, NUMROWS, devxout, devyout, devradprofileout, devimage, devbeadstatus, devradinprofile);
	
	cudaFree ( devradinprofile );
	
	//Prep the Radial Profiile and then fit it to the calibration stack to find the best fit height.
	PrepRadInProf(maxthreadsperblock, cross, NUMOFIMAGES, NUMOFBEADS, devradprofileout, devcosband, devCosBandLength);	

	//FitPrepandPhase<<<numBlocks,threadsPerBlock>>>(cross, devforgetradius, NUMOFIMAGES, NUMOFBEADS, devRealRow, devRealCol, devAmpRow, devAmpCol, devCorkRow, devCorkCol, devradprofileout, devReal, devamplitudes, devcorkscrews, devzout, devbeadstatus);
	FitPrepandPhase<<<blocksfitprep,threadsfitprep>>>(cross, devforgetradius, NUMOFIMAGES, NUMOFBEADS, devRealRow, devRealCol, devAmpRow, devAmpCol, devCorkRow, devCorkCol, devradprofileout, devReal, devamplitudes, devcorkscrews, devzout, devbeadstatus);
	
	//Copy the Z position back to Labview
	cudaMemcpy( (*zout)->val, devzout, NUMOFIMAGES*NUMOFBEADS*sizeof(float), cudaMemcpyDeviceToHost);	
//End Z Tracking

	//Check if the bead position for any of the beads are in a bad spot. If so, mark the bead as bad, and dont track it!
	CheckBead<<<numBlocks,threadsPerBlock>>>(NUMOFIMAGES, NUMOFBEADS, NUMROWS, NUMCOLS, cross, devRealRow, devbeadstatus, devxout, devyout, devzout);
	cudaMemcpy( (*BeadStatus)->val, devbeadstatus, NUMOFBEADS*sizeof(bool), cudaMemcpyDeviceToHost);

	//Free all allocated memory
	cudaFree ( devxfit ) ;
	cudaFree ( devyfit ) ;
	cudaFree ( devxdelta ) ;
	cudaFree ( devydelta ) ;
	cudaFree ( devAvgYProfileOut ) ;
	cudaFree ( devAvgXProfileOut ) ;
	cudaFree ( devAvgYProfileIn ) ;
	cudaFree ( devAvgXProfileIn ) ;
	cudaFree ( devcosband ) ;
	cudaFree ( devamplitudes );
	cudaFree ( devReal );
	cudaFree ( devcorkscrews );
	cudaFree ( devforgetradius );
	cudaFree ( devAmpRow );
	cudaFree ( devAmpCol );
	cudaFree ( devCorkRow );
	cudaFree ( devCorkCol );
	cudaFree ( devRealRow );
	cudaFree ( devRealCol );
	cudaFree ( devCosBandLength );
	cudaFree ( devzout ) ;
	cudaFree ( devradprofileout );
	cudaFree ( devimage ) ;
	cudaFree ( devyout ) ;
	cudaFree ( devxout ) ;
	cudaFree ( devyin ) ;
	cudaFree ( devxin ) ;
	cudaFree ( devbeadstatus ) ;

	//Check for any CUDA errors, and output it back to LabView.
	error = cudaGetLastError();
	if(error!=cudaSuccess) {
		memcpy(text, cudaGetErrorString(error), strlen(cudaGetErrorString(error)));
	} else
	{
		memcpy(text, "No_Errors", 9);
	}

	return 0;
}
