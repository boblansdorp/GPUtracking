//Using Complex numbers:
//complex in polar form is x = r * e(i* theta)
//abs(x) = r;
//arg(x) = theta;
void PhaseInNBHD (int cross, int forgetradius, int bestfit, double *phases, Array2dHandle amplitudes, complex<double> *radinprofile, CArray2dHandle corkscrews)
{
	int size = 5;
	int indexstart = (int)(bestfit - (size/2));
	int length = ((cross/2)-forgetradius);
	double temp;

	for(int i=0; i < size; i++)
	{
		double sum1=0, sum2=0;
		for(int j=0; j < length; j++)
		{
			temp = abs(radinprofile[j]) * (*amplitudes)->val[(((i+indexstart) * (*amplitudes)->length[1]) + j)];
			sum1 += temp;
			sum2 += temp * arg(radinprofile[j]/(*corkscrews)->val[(((i+indexstart) * (*corkscrews)->length[1]) + j)]);
		}
		phases[i]=sum2/sum1;
	}
}