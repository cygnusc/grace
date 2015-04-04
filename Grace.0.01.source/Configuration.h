#ifndef _CONFIGURATION_H_
#define _CONFIGURATION_H_

#include <string>
#include <iostream>
#include <fstream>
using namespace std;

class Configuration {
public:
	Configuration(std::string);
	~Configuration();
	bool readNext();
	int INDEX(int, int, int);

	ifstream infile;
	int frame;
	int totalFrames;

	float time;
	
	int dimx; // actual dimension on X direction
	int dimy;
	int dimz;

	int sampleX; // spacial sampling rate on X direction
	int sampleY;
	int sampleZ;
	
	int disDimX; // coordinates used in display arrays
	int disDimY;
	int disDimZ;

	int maxDimX;
	int maxDimY;
	int maxDimZ;

    int maxDim;

	float *X;
	float *Y;
	float *Z;

	float *vecX;
	float *vecY;
	float *vecZ;
};

#endif