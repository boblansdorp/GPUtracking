//This function corrects the XY output to correctly reflect its position on the image
__global__ void Out(int numofimages, int numofbeads, int cross, int *x, float *xout, bool *beadstatus)
{
	int ImageID = threadIdx.x + blockIdx.x * blockDim.x;
	int BeadID = threadIdx.y + blockIdx.y * blockDim.y;
if ((ImageID < numofimages) && (BeadID < numofbeads) && beadstatus[BeadID])
{
	xout[((ImageID*numofbeads) + BeadID)] -= (cross/2);
	xout[((ImageID*numofbeads) + BeadID)] /= 2;
	xout[((ImageID*numofbeads) + BeadID)] += x[BeadID];
} else if((ImageID < numofimages) && (BeadID < numofbeads) && !beadstatus[BeadID])
{
	xout[((ImageID*numofbeads) + BeadID)] = 0;
}
}