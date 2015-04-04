#ifndef UTIL_H
#define UTIL_H

#include <iostream>
#include <cassert>

class ExternH { // external magnetic field
public:
	float externHx;
	float externHy;
	float externHz;

	int startTime; // when the field appears
	int decayTime; // when the field begins to decay
	int endTime; // when the field vanishes

	static int numOfFields;

	ExternH(float, float, float, int, int, int);
	~ExternH();
	void print();
};



#endif
