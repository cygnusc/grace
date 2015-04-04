#include "util.h"

int ExternH::numOfFields = 0;

ExternH::ExternH(float ehx = 0.f, float ehy = 0.f, float ehz = 0.f, 
				 int sTime = 0, int dTime = 0,int eTime = 0) {
	externHx = ehx;
	externHy = ehy;
	externHz = ehz;
	startTime = sTime;
	decayTime = dTime;
	endTime = eTime;

	//numOfFields++;

	assert (startTime <= decayTime);
	assert (decayTime <= endTime);
}

ExternH::~ExternH() {
	//numOfFields--;
}

void ExternH::print() {
	std::cout << " external field = (" 
		<< externHx << " "
		<< externHy << " "
		<< externHz << "),"
		<< " time range = "
		<< startTime << " "
		<< decayTime << " "
		<< endTime << " "
		<< " Number of Fields = "
		<< numOfFields << " "
		<< std::endl;

}

ExternH *myExtH = new ExternH[2];