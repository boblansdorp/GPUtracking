//Defines 1d Array Handle
typedef struct {
	int length;
	double val[1];
		//to access value at a row => val[row]
} Array1d, **Array1dHandle;

//Defines 2d Array Handle
typedef struct {
	int length[2];
		//length[0] = number of rows
		//length[1] = number of columns
	double val[1];
		//to access value at a row/col => val[(row * length[1]) + col)]
} Array2d, **Array2dHandle;

//Defines 1d Complex Array Handle
typedef struct {
	int length;
	complex<double> val[1];
		//to access value at a row => val[row]
} CArray1d, **CArray1dHandle;

//Defines Upto 2d Complex Array Handle
typedef struct {
	int length[2]; //increase the index to handle bigger array sizes
		//length[0] = number of rows
		//length[1] = number of columns
	complex<double> val[1];
		//to access value at a row/col => val[(row * length[1]) + col)]
} CArray2d, **CArray2dHandle;

typedef struct{
	int x;
	int y;
	int cross;
	int crossthickness;
	Array2dHandle image;
} XYTrackInData;

typedef struct{
	//remember to make pointers
	double *xout;
	double *yout;
} XYTrackOutData;

typedef struct{
	int cross;
	int forgetradius;
	double newx;
	double newy;
	Array2dHandle image;
	Array2dHandle Real;
	Array2dHandle amplitudes;
	CArray1dHandle cosband;
	CArray2dHandle corkscrews;
} ZTrackInData;

typedef struct{
	int *index;
	double *z;
} ZTrackOutData;

typedef struct{
	int x;
	int y;
	int cross;
	int crossthickness;
	int forgetradius;
	Array2dHandle image;
	Array2dHandle Real;
	Array2dHandle amplitudes;
	CArray1dHandle cosband;
	CArray2dHandle corkscrews;
} XYZTrackInData;

typedef struct{
	int *index;
	double *xout;
	double *yout;
	double *z;
} XYZTrackOutData;