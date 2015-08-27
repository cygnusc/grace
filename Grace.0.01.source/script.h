/*

script.h - Grace Micromagnetic Simulator written by Ru Zhu, zhu@sting.graceland.edu

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

#ifndef _SCRIPT_H_
#define _SCRIPT_H_

#include <fstream>
#include <string>

class Script {
public:
    int readScript(int argc, char** argv);
    Script();
    void print();
    void print(std::ofstream&);
    void readRunID(std::fstream &);
    void readGlobalField(std::fstream &);
    void readExternField(std::fstream &);
    void readRectang(std::fstream &);
    void readMaterial(std::fstream &);
    void readSimulation(std::fstream &);
    void readOutput(std::fstream &);

    std::string runID;

    int nx;
    int ny;
    int nz;

    float Ms;
    float alpha;
    float exchConstant;
    float initMx;
    float initMy;
    float initMz;

    float Hkx;
    float Hky;
    float Hkz;

	bool readInitState;
	bool writeFinalState;

    int writeInterval;
    int timesteps;
    float dt;

    float externFieldx;
    float externFieldy;
    float externFieldz;

	float externFieldStartTime;
	float externFieldDecayTime;
	float externFieldStopTime;

};

#endif
