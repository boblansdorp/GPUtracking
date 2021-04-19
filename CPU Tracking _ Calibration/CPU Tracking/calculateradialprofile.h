//This program will analyze 'all' pixels of our subimage, and create a radial profile of the intensity
//'all' is in quotes because we actually only care about pixels that will form full and complete circles thus corners will be ignored to a certain extent

//NOTE we have a bit of discrepency in the last values for the radial profile
//also we might want to include a check for acessing image data out of bounds.
void CalculateRadialProfile(int cross, double x, double y, double *radprofileout, Array2dHandle image)
{
	int numcol = (*image)->length[1];

	//This is where our radial profile will be stored
	//radinprofile[radius][0] = sum of intensity
	//radinprofile[radius][1] = sum of partial weight
	double **radinprofile;
	radinprofile = new double*[(cross/2)];
	for(int i = 0; i < (cross/2); i++)
	{
		radinprofile[i] = new double[2];
		radinprofile[i][0] = 0;
		radinprofile[i][1] = 0;
	}

	//we start analyzing at cross/2 before the central x and y position
	//ceil rounds up to next highest integer
	int xstart = (int)(ceil(x) - (cross/2));
	int ystart = (int)(ceil(y) - (cross/2));

	double xfrac = x - ceil(x);
	double yfrac = y - ceil(y);

	for(int i = 0; i < cross; i++)
	{
		for(int j = 0; j < cross; j++)
		{
			//note that we must subtract the center location, cross/2 from the index to calulate radius from the central location
			double radius = sqrt((pow(((cross/2)+yfrac-i),2) + pow(((cross/2)+xfrac-j),2)));
			//skip points on the corners that do not form complete circles
			if (radius < (cross/2))
			{
				double radfloor;
				double radfrac = modf(radius, &radfloor);
				double radceil = radfloor + 1;
				//Add a weighted intensity profile and the weighting value
				radinprofile[(int)radfloor][0] += ((1-radfrac) * (*image)->val[((ystart+i) * numcol) + (xstart+j)]);
				radinprofile[(int)radfloor][1] += (1-radfrac);
				//this fixes an error when radius is 63.X and then it tried to assign a value to radinprofile[64][0] which doesnt exist
				if (radceil < (cross/2))
				{
					radinprofile[(int)radceil][0] += (radfrac * (*image)->val[((ystart+i) * numcol) + (xstart+j)]);
					radinprofile[(int)radceil][1] += radfrac;
				}
			}
		}
	}
	
	for(int i = 0; i < (cross/2); i++)
	{
		//divide the intensity profile by the total weight to normalize the profile
		if (radinprofile[i][1] > 0)
			radprofileout[i] = (radinprofile[i][0]/radinprofile[i][1]);
		else
			radprofileout[i] = 0;
		delete[] radinprofile[i];
	}
	delete[] radinprofile;
}