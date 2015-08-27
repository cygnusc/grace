/*

Configuration.h - Grace Micromagnetic Simulator written by Ru Zhu, zhu@sting.graceland.edu

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

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
