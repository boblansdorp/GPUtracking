//square vector with cosine bandpass
__global__ void SquareVector( int cross, int numofimages, int numofbeads, cufftComplex *a)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	while (tid < cross * numofimages * numofbeads)
	{
		if(tid%cross < 60)
		{
			a[tid].x *= (1+cos((tid%cross)*(2*PI/cross)))/2;
			a[tid].y *= (1+cos((tid%cross)*(2*PI/cross)))/2;
		} else
		{
			a[tid].x *= 0;
			a[tid].y *= 0;
		}
		float tempREAL = a[tid].x;
		a[tid].x = a[tid].x * a[tid].x - (a[tid].y *a[tid].y);
		a[tid].y = 2*tempREAL*a[tid].y;
		tid += blockDim.x * gridDim.x;
	}
}