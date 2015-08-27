/*

util.h - Grace Micromagnetic Simulator written by Ru Zhu, zhu@sting.graceland.edu

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
