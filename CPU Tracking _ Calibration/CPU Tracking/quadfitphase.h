//This takes a descrete value for the best fit cal image index, and then uses a quadratic fit
//ong the phase information to find a non integer more accurate value for the best fit z height

double QuadFitPhase (int bestfit, double *phases)
{
	int size = 5;
	double **data;
	double solution[3] = {0};
	data = new double*[size];
	for(int i=0; i < size; i++)
	{
		data[i] = new double[2];
		//note that phases is the domain
		data[i][0] = phases[i];
		data[i][1] = (i - floor((double)size/2));
	}

	QuadFit(size, data, solution);
	return (bestfit + solution[2]);
	for (int i=0; i < size; i++)
		delete[] data[i];
	delete[] data;
}