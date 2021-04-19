//Yay linear Algebra! This is my attempt at solving a quadratic function using 5 data points.
//General theory: With points [x1,y1]...[x5,y5], and trying to fit to y = ax^2 + bx + c
//In matrix form this becomes: a = [a,b,c]^T, y = [y1,...,y5]^T, and X = [(x1^2, x1, 1),...,(x5^2, x5, 1)]
//In General... => (X)a = y => (X^T * X)a = (X^T)y => a = ((X^T * X)^(-1) * (X^T))y
//(X^T * X) = [(sum(xi^4), sum(xi^3), sum(xi^2)), (sum(xi^3), sum(xi^2), sum(xi)), (sum(xi^2), sum(xi), sum(i))]
//det((X^T *X)) = sum(xi^4)(sum(xi^2)sum(i)-sum(xi)^2) + sum(xi^3)(sum(xi)sum(xi^2)-sum(i)sum(xi^3)) + sum(xi^2)(sum(xi^3)sum(xi)-sum(xi^2)^2)
// => (X^T *X)^-1 = 1/det((X^T *X)) * [{sum(i)sum(xi^2) - sum(xi)^2, -(sum(i)sum(xi^3) - sum(xi)sum(xi^2)), sum(xi)sum(xi^3)-sum(xi^2)^2}, {-(sum(i)sum(xi^3)-sum(xi^2)sum(xi), sum(i)sum(xi^4)-sum(xi^2)^2, -(sum(xi)sum(xi^4)-sum(xi^3)sum(xi^2))}, {sum(xi)sum(xi^3) - sum(xi^2)sum(xi^2), -(sum(xi)sum(xi^4) - sum(xi^2)sum(xi^3)), sum(xi^2)sum(xi^4) - sum(xi^3)sum(xi^3)}]
// => (X^T *X)^-1 *X^T = stuff...
//Using this info, we can now solve for coefficents!

void QuadFit(int size, double **data, double *a)
{
	//size determines how many points used
	//data should be a 2xsize matrix with x value in first place and y value in second.
	//data[3][0] -> x4, data[3][1] -> y4

	// sumx3 = (x0)^3 + (x1)^3 + ... + (x{size})^3
	double sumx4 = 0, sumx3 = 0, sumx2 = 0, sumx1 = 0, det;
	
	// Allocate memory for X^T
	double **transpose;
	transpose = new double*[3];
	for(int i=0; i < 3; i++)
		transpose[i] = new double[size];

	//Allocate memory for (X^T *X)^-1 * X^T
	double **stuff;
	stuff = new double*[3];
	for(int i=0; i < 3; i++)
		stuff[i] = new double[size];

	//input data, do some precalculations at the same time to avoid making more loops
	for (int i=0; i < size; i++)
	{
			//these computations only happen with x values
			int j=0;
			{
				//develops the sum values needed for inverse
				sumx4 += data[i][j]*data[i][j]*data[i][j]*data[i][j];
				sumx3 += data[i][j]*data[i][j]*data[i][j];
				sumx2 += data[i][j]*data[i][j];
				sumx1 += data[i][j];
				
				//develops transpose matrix
				transpose[2][i] = 1;
				transpose[1][i] = data[i][j];
				transpose[0][i] = data[i][j]*data[i][j];
			}
	}
	
//After solving all the math
	//determinate
	det = (sumx4*sumx2*size) + (sumx3*sumx1*sumx2) + (sumx2*sumx3*sumx1) - (sumx2*sumx2*sumx2) - (sumx1*sumx1*sumx4) - (size*sumx3*sumx3); 

	//precalculated the inverse matrix to avoid numerical methods which take time or lose accuracy, NOTE: does not include division of determinate
	double inverse[3][3] = {
		{size*sumx2 - sumx1*sumx1, -(size*sumx3 - sumx1*sumx2), sumx1*sumx3-sumx2*sumx2},
		{-(size*sumx3-sumx2*sumx1), size*sumx4-sumx2*sumx2, -(sumx1*sumx4-sumx3*sumx2)},
		{sumx1*sumx3 - sumx2*sumx2, -(sumx1*sumx4 - sumx2*sumx3), sumx2*sumx4 - sumx3*sumx3}
	};

	//This is matrix multiplication for this particular pair of matrices
	for (int i=0; i < 3; i++)
	{
		for (int j=0; j < size; j++)
		{
			stuff[i][j] = inverse[i][0]*transpose[0][j] + inverse[i][1]*transpose[1][j] + inverse[i][2]*transpose[2][j];
		}
	}

	//This is the final matrix multiplication that outputs a 1x3 matrix with our curve parameters
	for (int i=0; i < 3; i++)
	{
		for (int j=0; j < size; j++)
		{
			a[i] += stuff[i][j]*data[j][1];
		}
		//dont forget to divide by determinate
		a[i] /= det;
	}
}