#define PI 3.14159265
#define DllExport __declspec(dllexport)

//We must use this so that data is packaged correctly for LabView
//In 64 bit Labview apparently we don't need this
//#pragma pack(1)

//Only Necessary for the standalone EXE File
#include <iostream>
#include <fstream>
using namespace std;
//Necessary Includes for Math
#include <cuda.h>
#include <cufft.h>
#include <complex>
#include <math.h>
#include <algorithm>

//Other
#include "checkbead.h"
#include "quadfit.h"
//XY Tracking Header Files
#include "datatypes.h"
#include "avgycrossandprep.h"
#include "avgxcrossandprep.h"
#include "squarevector.h"
#include "finddelta.h"
#include "rotate.h"
#include "out.h"
#include "removediscretization.h"
//Z Tracking Header Files
#include "calculateradialprofile1.h"
#include "calculateradialprofile2.h"
#include "prepradprofile.h"
#include "fitprepandphase.h"
//Main Tracking Program
#include "GPUTracking.h"