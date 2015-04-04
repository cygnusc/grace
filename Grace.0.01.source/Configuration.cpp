#include "Configuration.h"
#include <cassert>
#define FILE_OPEN_ERROR -1
#define READ_DONE 0
#define TEST 0
//#define INDEX(i,j,k) (k*dimx*dimy + j*dimx+ i)

int Configuration::INDEX(int i, int j, int k) {
    return (k*disDimX*disDimY + j*disDimX+ i);
}

Configuration::Configuration(string filename) {
	infile.open(filename);
	if (!infile.good()) {
		cout << " cannot open file " << filename << ", press ENTER to exit" << endl;
		getchar(); 
		exit(FILE_OPEN_ERROR);
	}
    char junkline [1024];
	while (infile.peek() == '#') {
        infile.getline(junkline, 1024); // get rid of file headers
	}

	string junk;
	infile >> junk;
	assert (junk == "dims");
	infile >> dimx >> dimy >> dimz; // get dimensions
	cout << "dimensions = " << dimx << " " << dimy << " " << dimz << endl;

	

	infile >> junk >> sampleX >> sampleY >> sampleZ;
	assert (junk == "sampling");
	cout << "sampling = " << sampleX << " " << sampleY << " " << sampleZ << endl;

	disDimX = dimx / sampleX;
	disDimY = dimy / sampleY;
	disDimZ = dimz / sampleZ;
	cout << "disDims = " << disDimX << " " << disDimY << " " << disDimZ << " " << endl;

	X = new float [disDimX * disDimY * disDimZ];
	Y = new float [disDimX * disDimY * disDimZ];
	Z = new float [disDimX * disDimY * disDimZ];

	vecX = new float [disDimX * disDimY * disDimZ];
	vecY = new float [disDimX * disDimY * disDimZ];
	vecZ = new float [disDimX * disDimY * disDimZ];

	infile >> junk >> totalFrames;
	frame = totalFrames;
	assert (junk == "frames");
	cout << "frame = " << frame << endl;

	maxDimX = dimx - sampleX;
	maxDimY = dimy - sampleY;
	maxDimZ = dimz - sampleZ;

    maxDim = 0;
    maxDim = maxDimX > maxDimY? maxDimX : maxDimY;
    maxDim = maxDim > maxDimZ ? maxDim : maxDimZ;
#if TEST
	cout << " max Dim X, Y, Z = " << maxDimX << " " << maxDimY << " " << maxDimZ << endl;
#endif
}

Configuration::~Configuration() {
	infile.close();
}

bool Configuration::readNext() {
	float trash;
	if (infile.good() == false) {
		cout << " cannot read file, press ENTER to exit" << endl;
		getchar();
		exit(FILE_OPEN_ERROR);
	}
	if (frame == 0) {
		cout << " reaches the end of the file." << endl;
		exit (READ_DONE);
	}
	char junk;
	infile >> junk >> time;
	assert(junk == 't');
#if TEST
	cout << "time = " << time << endl;
#endif
	for (int i = 0; i < disDimX; i+=1) {
	for (int j = 0; j < disDimY; j+=1) {
	for (int k = 0; k < disDimZ; k+=1) {
#if TEST
		cout << k << " " << j << " " << i << " ";
#endif
		infile >> X[INDEX(i,j,k)] >> Y[INDEX(i,j,k)] >> Z[INDEX(i,j,k)];
//		if (X[INDEX(i,j,k)] > maxDimX) maxDimX = X[INDEX(i,j,k)]; 
//		if (Y[INDEX(i,j,k)] > maxDimY) maxDimY = Y[INDEX(i,j,k)]; 
//		if (Z[INDEX(i,j,k)] > maxDimZ) maxDimZ = Z[INDEX(i,j,k)]; 
#if TEST
		cout << X[INDEX(i,j,k)] << " " << Y[INDEX(i,j,k)] << " " << Z[INDEX(i,j,k)] << " ";
#endif
		infile >> vecX[INDEX(i,j,k)] >> vecY[INDEX(i,j,k)] >> vecZ[INDEX(i,j,k)] >> trash >> trash >> trash;
#if TEST
		cout << vecX[INDEX(i,j,k)] << " " << vecY[INDEX(i,j,k)] << " " << vecZ[INDEX(i,j,k)] << endl;
#endif
	}
	}
	}


	char garbage[256];
	infile.getline(garbage, 256);
	frame --;
	return true;
}

Configuration * myConfig; //("example.txt");