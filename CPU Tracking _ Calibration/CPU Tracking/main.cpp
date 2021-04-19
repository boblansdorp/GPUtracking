//Use this to fix the way double packets are sent to and from labview
//#pragma pack(1)

//Using standard names and such, for some reason I MUST include iostream although I dont use it.
#include <iostream>
using namespace std;

//Used for Rotation of array
#include <algorithm>

//Allows the use of complex numbers and such
//Be sure to include this BEFORE FFTW3 so that fftw_complex just becomes a normal complex number
#include <complex>

//Fast Fourier Transform Header File
#include "fftw3.h"
//Allows Labview to send in and recieve data
#define DllExport __declspec(dllexport)

//Broke up the program into header files, these are the components
#include "datatypes.h"
#include "findavgprofilecenter.h"
#include "prepavgprofile.h"
//#include "finddelta.h" is a part of findavgprofilecenter
////#include "quadfit.h" is a part of finddelta
#include "avgxycross.h"

//These are for z tracking
#include "calculateradialprofile.h"
#include "quadfitphase.h"
#include "fitpreptocal.h"
#include "phaseinnbhd.h"
#include "prepradinprof.h"

//main program
//DllExport void CPUTracking (XYZTrackInData globalin, XYZTrackOutData globalout)
DllExport void CPUTracking (int x, int y, int cross, int crossthickness, int forgetradius, Array2dHandle image, Array2dHandle Real, Array2dHandle amplitudes, CArray1dHandle cosband, CArray2dHandle corkscrews, int *index, double *xout,	double *yout, double *z)
{
	//XYZTrackInData xyztrackin = globalin;
	//int out1 = 0;
	//double out2=0, out3=0, out4=0;
	//XYZTrackOutData xyztrackout = {&out1,&out2,&out3,&out4};

	//These variables are passed between functions, but never to labview.
	//double *radinprofile = new double[(xyztrackin.cross/2)];
	double *radinprofile = new double[(cross/2)];

	//double *prep = new double[((xyztrackin.cross/2)-xyztrackin.forgetradius)];
	double *prep = new double[((cross/2)-forgetradius)];
	//complex<double> *corkscrew = new complex<double>[((xyztrackin.cross/2)-xyztrackin.forgetradius)];
	complex<double> *corkscrew = new complex<double>[((cross/2)-forgetradius)];
	double *phases = new double[5];

	//double *AvgXProfile = new double[xyztrackin.cross];
	double *AvgXProfile = new double[cross];
	//double *AvgYProfile = new double[xyztrackin.cross];
	double *AvgYProfile = new double[cross];

	//Start XY Tracking
	//AvgXYCross(xyztrackin.x, xyztrackin.y, xyztrackin.cross, xyztrackin.crossthickness, xyztrackin.image, AvgXProfile, AvgYProfile);
	AvgXYCross(x, y, cross, crossthickness, image, AvgXProfile, AvgYProfile);
	//PrepAvgProfile (xyztrackin.cross, AvgXProfile, AvgYProfile);
	
	PrepAvgProfile (cross, AvgXProfile, AvgYProfile);
	
//	FindAvgProfileCenter (xyztrackin.cross, xyztrackin.x, xyztrackin.y, AvgXProfile, AvgYProfile, xyztrackout.xout, xyztrackout.yout);
	
	FindAvgProfileCenter (cross, x, y, AvgXProfile, AvgYProfile, xout, yout);

	//End XY Tracking
	
	//Start Z Tracking
	//CalculateRadialProfile(xyztrackin.cross, *xyztrackout.xout, *xyztrackout.yout, radinprofile, xyztrackin.image);
	CalculateRadialProfile(cross, *xout, *yout, radinprofile, image);
	//PrepRadInProf (xyztrackin.cross, xyztrackin.forgetradius, radinprofile, prep, corkscrew, xyztrackin.cosband);
	PrepRadInProf (cross, forgetradius, radinprofile, prep, corkscrew, cosband);
	//*xyztrackout.index = FitPreptoCal (prep, xyztrackin.Real);
	*index = FitPreptoCal (prep, Real);
	//PhaseInNBHD (xyztrackin.cross, xyztrackin.forgetradius, *xyztrackout.index, phases, xyztrackin.amplitudes, corkscrew, xyztrackin.corkscrews);
	PhaseInNBHD (cross, forgetradius, *index, phases, amplitudes, corkscrew, corkscrews);
	//*xyztrackout.z = QuadFitPhase (*xyztrackout.index, phases);
	*z = QuadFitPhase (*index, phases);
	//End Z Tracking

	//*globalout.index = *xyztrackout.index;
	//*globalout.xout = *xyztrackout.xout;
	//*globalout.yout = *xyztrackout.yout;
	//*globalout.z = *xyztrackout.z;


	//Free Memory
	delete[] radinprofile;
	delete[] prep;
	delete[] corkscrew;
	delete[] phases;
	delete[] AvgXProfile;
	delete[] AvgYProfile;
	
}

DllExport void Calibration (int x, int y, int cross, int crossthickness, Array2dHandle image, double *xout, double *yout, Array1dHandle RadProfileOut)
{
	double *AvgXProfile = new double[cross];
	double *AvgYProfile = new double[cross];
	double *radinprofile = new double[(cross/2)];
	//Start XY Tracking
	AvgXYCross(x, y, cross, crossthickness, image, AvgXProfile, AvgYProfile);
	PrepAvgProfile (cross, AvgXProfile, AvgYProfile);
	FindAvgProfileCenter (cross, x, y, AvgXProfile, AvgYProfile, xout, yout);

	//End XY Tracking
	CalculateRadialProfile(cross, *xout, *yout, radinprofile, image);
	
	for(int i=0; i < cross/2; i++)
		(*RadProfileOut)->val[i] = radinprofile[i];

	delete[] AvgXProfile;
	delete[] AvgYProfile;
	delete[] radinprofile;

}

DllExport void RadialProfile (int cross, float x, float y, Array1dHandle RadProfileOut, Array2dHandle image)
{
	double *radinprofile = new double[(cross/2)];

	//End XY Tracking
	CalculateRadialProfile(cross, x, y, radinprofile, image);
	
	for(int i=0; i < cross/2; i++)
		(*RadProfileOut)->val[i] = radinprofile[i];

	delete[] radinprofile;
}