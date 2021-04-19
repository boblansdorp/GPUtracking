#include "include.h"

//This Whole Function allows us to run the tracking program independent of LabView.
//Use our Labview Program to Output the tracking data into a file, and the modify the filename below.
//To modify the tracking program look at GPUTracking.h
int main()
{
	Array1dIntHandle x = new Array1dInt*;
	Array1dIntHandle y = new Array1dInt*;
	int cross;
	int crossthickness;
	ImageHandle image = new Image*;
	Array2dHandle xout = new Array2d*;
	Array2dHandle yout = new Array2d*;
	Array2dHandle zout = new Array2d*;
	ArrayClusterHandle calibration = new ArrayCluster*;
	BoolArrayHandle BeadStatus = new BoolArray*;
	char* text = new char[1000];
	Array1dHandle test = new Array1d*;

	ifstream file ("TESTDATAINNOTEXTforget.dat", ios::in|ios::binary);

	if(file.is_open())
	{	
		int tempsize[3] = {0};
		file.read((char *)&tempsize[0], sizeof(int));
		(*x) = (Array1dInt*)malloc(sizeof(int) + sizeof(int)*tempsize[0]);
		(*x)->length = tempsize[0];
		file.read((char *)&(*x)->val[0], sizeof(int)*tempsize[0]);

		file.read((char *)&tempsize[0], sizeof(int));
		(*y) = (Array1dInt*)malloc(sizeof(int) + sizeof(int)*tempsize[0]);
		(*y)->length = tempsize[0];
		file.read((char *)&(*y)->val[0], sizeof(int)*tempsize[0]);

		file.read((char *)&cross, sizeof(int));
		file.read((char *)&crossthickness, sizeof(int));

		file.read((char *)&tempsize[0], sizeof(int)*3);
		(*image) = (Image*)malloc(sizeof(int)* 3 + sizeof(unsigned char)*tempsize[0]*tempsize[1]*tempsize[2]);
		for(int i=0; i < 3; i++)
			(*image)->length[i] = tempsize[i];
		file.read((char *)&(*image)->val[0], sizeof(unsigned char)*tempsize[0]*tempsize[1]*tempsize[2]);

		file.read((char *)&tempsize[0], sizeof(int)*2);
		(*xout) = (Array2d*)malloc(sizeof(int)* 2 + sizeof(float)*tempsize[0]*tempsize[1]);
		for(int i=0; i < 2; i++)
			(*xout)->length[i] = tempsize[i];
		file.read((char *)&(*xout)->val[0], sizeof(float)*tempsize[0]*tempsize[1]);

		file.read((char *)&tempsize[0], sizeof(int)*2);
		(*yout) = (Array2d*)malloc(sizeof(int)* 2 + sizeof(float)*tempsize[0]*tempsize[1]);
		for(int i=0; i < 2; i++)
			(*yout)->length[i] = tempsize[i];
		file.read((char *)&(*yout)->val[0], sizeof(float)*tempsize[0]*tempsize[1]);

		file.read((char *)&tempsize[0], sizeof(int)*2);
		(*zout) = (Array2d*)malloc(sizeof(int)* 2 + sizeof(float)*tempsize[0]*tempsize[1]);
		for(int i=0; i < 2; i++)
			(*zout)->length[i] = tempsize[i];
		file.read((char *)&(*zout)->val[0], sizeof(float)*tempsize[0]*tempsize[1]);
		
		//CLUSTER START
		file.read((char *)&tempsize[0], sizeof(int));
		(*calibration) = (ArrayCluster*)malloc(sizeof(int) + sizeof(Cluster) * tempsize[0]);
		(*calibration)->length = tempsize[0];
		for (int i=0; i< tempsize[0]; i++)
		{	
			int temp[2] = {0};
			//forgetradius
			file.read((char *)&(*calibration)->cluster[i].forgetradius, sizeof(int));
			//zstep
			file.read((char *)&(*calibration)->cluster[i].zstep, sizeof(float));
			//cosband
			(*calibration)->cluster[i].cosband = new Array1d*;
			file.read((char *)&temp[0], sizeof(int));
			(*(*calibration)->cluster[i].cosband) = (Array1d*)malloc(sizeof(int) + sizeof(float) * temp[0]);
			(*(*calibration)->cluster[i].cosband)->length = temp[0];
			file.read((char *)&(*(*calibration)->cluster[i].cosband)->val[0], sizeof(float) * temp[0]);
			//Amplitude
			(*calibration)->cluster[i].amplitude = new Array2d*;
			file.read((char *)&temp[0], sizeof(int)*2);
			(*(*calibration)->cluster[i].amplitude) = (Array2d*)malloc(sizeof(int)*2 + sizeof(float)*temp[0]*temp[1]);
			for(int j=0; j < 2; j++)
				(*(*calibration)->cluster[i].amplitude)->length[j] = temp[j];
			file.read((char *)&(*(*calibration)->cluster[i].amplitude)->val[0], sizeof(float) * temp[0]*temp[1]);
			//corkscrews
			(*calibration)->cluster[i].corkscrews = new CArray2d*;
			file.read((char *)&temp[0], sizeof(int)*2);
			(*(*calibration)->cluster[i].corkscrews) = (CArray2d*)malloc(sizeof(int)*2 + sizeof(complex<float>)*temp[0]*temp[1]);
			for(int j=0; j < 2; j++)
				(*(*calibration)->cluster[i].corkscrews)->length[j] = temp[j];
			file.read((char *)&(*(*calibration)->cluster[i].corkscrews)->val[0], sizeof(complex<float>) * temp[0]*temp[1]);
			//Real
			(*calibration)->cluster[i].Real = new Array2d*;
			file.read((char *)&temp[0], sizeof(int)*2);
			(*(*calibration)->cluster[i].Real) = (Array2d*)malloc(sizeof(int)*2 + sizeof(float)*temp[0]*temp[1]);
			for(int j=0; j < 2; j++)
				(*(*calibration)->cluster[i].Real)->length[j] = temp[j];
			file.read((char *)&(*(*calibration)->cluster[i].Real)->val[0], sizeof(float) * temp[0]*temp[1]);
		}
		//CLUSTER END

		file.read((char *)&tempsize[0], sizeof(int));
		(*BeadStatus) = (BoolArray*)malloc(sizeof(int) + sizeof(bool)*tempsize[0]);
		(*BeadStatus)->length = tempsize[0];
		file.read((char *)&(*BeadStatus)->val[0], sizeof(bool)*tempsize[0]);

		file.read((char *)&tempsize[0], sizeof(int));
		cout << tempsize[0] << endl;
		(*test) = (Array1d*)malloc(sizeof(int) + sizeof(float)*tempsize[0]);
		(*test)->length = tempsize[0];
		file.read((char *)&(*test)->val[0], sizeof(float)*tempsize[0]);
		
		cout << file.tellg() << endl;
		file.seekg(0, ios_base::end);
		cout << file.tellg() << endl << endl;

		GPUTracking (0, 1, 1, 1, 1, x, y, cross, crossthickness, image, xout, yout, zout, calibration, BeadStatus,text,test);

		for (int i=0; i < (*xout)->length[0] * (*xout)->length[1]; i++)
		{
			cout << i/2;
			if (i%2 == 0)
				cout << "a";
			else
				cout << "b";
			cout << ": "<< (*xout)->val[i] << " | " << (*yout)->val[i] << " | " << (*zout)->val[i] << endl;
		}
		cout << text << endl;
		
		
		system("pause");
		//free them all!
		free((*x));
		free((*y));
		free((*image));
		free((*xout));
		free((*yout));
		free((*zout));
		free((*calibration));
		free((*BeadStatus));
		free((*test));
		file.close();
	} else
	{
		cout << "File Did Not Open!" << endl;
	}

	delete x;
	delete y;
	delete image;
	delete xout;
	delete yout;
	delete zout;
	delete calibration;
	delete BeadStatus;
	delete text;
	delete test;

	return 0;
}