/*

util.cpp - Grace Micromagnetic Simulator written by Ru Zhu, zhu@sting.graceland.edu

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
