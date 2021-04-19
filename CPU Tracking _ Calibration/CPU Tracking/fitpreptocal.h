//this finds the least squares fit of the prepped I(r) to find the index of the closest slice
int FitPreptoCal (double *Prep, Array2dHandle Real)
{
	int index=0;
	int NumRow = (*Real)->length[0];
	int NumCol = (*Real)->length[1];
	double *Arr = new double[NumRow];
	for (int i=0; i < NumRow; i++)
	{
		Arr[i]=0;
		for (int j=0; j < NumCol; j++)
		{
			Arr[i] += pow((Prep[j] - (*Real)->val[((i * NumCol) + j)]), 2);
		}
		if(Arr[i] < Arr[index])
			index = i;
	}
	return index;
	delete[] Arr;
}