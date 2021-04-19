//This Preps the Data so that it can be Fourier Transformed
void PrepAvgProfile (int cross, double *AvgXProfile, double *AvgYProfile)
{
	//This brings the graph down to zero by finding the average height under the first and last quarter of the graph, 
	//and then subtracting all points by this amount.
	double Xelvfix = 0;
	double Yelvfix =0;
	//Adds all the elevation from the first and last quarter
	for(int i = 0; i < cross/4; i++)
	{
		Xelvfix += AvgXProfile[i] + AvgXProfile[cross - 1 - i];
		Yelvfix += AvgYProfile[i] + AvgYProfile[cross - 1 - i];
	}
	//Finds Average Elevation
	Xelvfix /= (cross/2);
	Yelvfix /= (cross/2);
	//Subtracks all values by the average elevation
	for(int i = 0; i < cross; i++)
	{
		AvgXProfile[i] = (AvgXProfile[i] - Xelvfix) * (1-cos(2*(3.14159)*i/(cross-1)))/2;
		AvgYProfile[i] = (AvgYProfile[i] - Yelvfix) * (1-cos(2*(3.14159)*i/(cross-1)))/2;
	}
}