/*

display.h - Grace Micromagnetic Simulator written by Ru Zhu, zhu@sting.graceland.edu

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

#ifndef _DISPLAY_H_
#define _DISPLAY_H_

#define DISPLAY_TEST 0
#define PI 3.14159265f

class Display{
public:

	void Run(int, char**);

};

void GLInit(int, char**);
void idleDisplay(void);
void reShape(int, int);
void displayVec(float, float, float, float, float, float);
void displayFunc();
void keyboardFunc(unsigned char, int, int);
void mouseButton(int, int, int, int);
void mouseMove(int, int);


#endif
