//this will make a profile of the bead using a + cross along the x and y axis at the approximate center of the bead

void AvgXYCross(int x, int y, int cross, int crossthickness, Array2dHandle image, double *AvgXProfile, double *AvgYProfile)
{
	//Initial Declarations
	int numcol = (*image)->length[1];

	//Starting from the center, find the top left corner of subarray.
	int RowStart = (int)(y - (cross * .5));
	int ColStart = (int)(x - (cross * .5));
	int YSliceStart = (int)(y - (crossthickness * .5));
	int XSliceStart = (int)(x - (crossthickness * .5));
	//Sum Slices Along Y -> AvgYProfile
	for(int j = 0; j < cross; j++)
	{
		AvgXProfile[j] = 0;
		for (int i = 0; i < crossthickness; i++)
		{
			AvgXProfile[j] += (*image)->val[((i+YSliceStart) * numcol) + (j+ColStart)];
		}
		//Finds the Average
		AvgXProfile[j] /= crossthickness;
	}
	//Sum Slices Along X -> AvgXProfile
	//-Notice i and j flipped in loops, but i still means Rows and j means Cols
	for(int i = 0; i < cross; i++)
	{
		AvgYProfile[i] = 0;
		for (int j = 0; j <crossthickness; j++)
		{
			AvgYProfile[i] += (*image)->val[((i+RowStart) * numcol) + (j+XSliceStart)];
		}
		AvgYProfile[i] /= crossthickness;
	}
}