#include "script.h"
#include <cstring>
#include <iostream>
#include <cassert>
#include "util.h"
//#include <fstream>

extern ExternH *myExtH;

int Script::readScript(int argc, char **argv) {
	std::fstream scriptFile;
	char* filename = new char[256];
	char* instr = new char[1024];
	char* junk = new char[1024];
	//char c;
	if (argc == 1) { // no argument provided
		filename = "default.txt";
	} else { // user defined script file provided
		filename = argv[1];
	}
	scriptFile.open(filename);

	runID = filename;

	runID = runID.substr(0, runID.length() - 4);
	std::cout << "runID = " << runID << std::endl;

	if (!scriptFile.is_open()) {
		std::cout << "Can't open script file " << filename
		<< ". Please check if it exists or being used by another program."
		<< std::endl;
		exit(1);
	}
	do {
		
			scriptFile >> instr;
			//std::cout << instr << std::endl;
            
			/*if (!strcmp(instr, "-runID")) {
                std::cout << "Reading runID info... " << std::endl;
                readRunID(scriptFile);
			} else 
			*/
			if (!strcmp(instr, "-globalfield")) {
				std::cout << "Reading global field... " << std::endl;
				readGlobalField(scriptFile);
			} else if (!strcmp(instr, "-externfield")) {
				std::cout << "Reading extern field... " << std::endl;
				readExternField(scriptFile);
			} else if (!strcmp(instr, "-rectang")) {
				std::cout << "Reading rectangular info... " << std::endl;
				readRectang(scriptFile);
			} else if (!strcmp(instr, "-readInitState")) {
				std::cout << "ReadInitState set to true... " << std::endl;
				readInitState = true;
			} else if (!strcmp(instr, "-writeFinalState")) {
				std::cout << "will write final state to " << runID << ".fc.txt" << std::endl;
				writeFinalState = true;
			} else if (!strcmp(instr, "-material")) {
				std::cout << "Reading material info... " << std::endl;
				readMaterial(scriptFile);
			} else if (!strcmp(instr, "-simulation")) {
				std::cout << "Reading simulation info... " << std::endl;
				readSimulation(scriptFile);
			} else if (!strcmp(instr, "-output")) {
				std::cout << "Reading output info... " << std::endl;
				readOutput(scriptFile);
			} else if (instr[0] == '#' ) {
				//std::cout << "skipping comment line... " << std::endl;
				scriptFile.getline(junk, 1024);
				//std::cout << "junk skipped = " << junk << std::endl;
			} else if (instr[0] != '\0') {
				std::cout << "Line not understood... " << std::endl;
				std::cout << instr << std::endl;
			}
	}
	while (!scriptFile.eof());

	scriptFile.close();

	return 0;
}

Script::Script() {
    runID = "";
	nx = ny = nz = 16;

	Ms = 1000.f;
	alpha = 1.f;
	exchConstant = 1.E-11f * 1.E18f;
	initMx = 1000.f;
	initMy = 0.f;
	initMz = 0.f;

	Hkx = 0.f;
	Hky = 0.f;
	Hkz = 0.f;

	readInitState = false;
	writeFinalState = false;
	
	writeInterval = 1;
	timesteps = 100;
	dt = 0.00001f;

	externFieldx = 0.f;
	externFieldy = 0.f;
	externFieldz = 0.f;
}

void Script::print() {
	std::cout << "nx = " << nx << std::endl
          << "ny = " << ny << std::endl
          << "nz = " << nz << std::endl
		  << "Ms = " << Ms << std::endl
		  << "alpha = " << alpha << std::endl
		  << "exchConstant = " << exchConstant << std::endl
		  << "initMx = " << initMx << std::endl
		  << "initMy = " << initMy << std::endl
		  << "initMz = " << initMz << std::endl
		  << "Hkx = " << Hkx << std::endl
		  << "Hky = " << Hky << std::endl
		  << "Hkz = " << Hkz << std::endl
		  << "writeInterval = " << writeInterval << std::endl
		  << "timesteps = " << timesteps << std::endl
		  << "dt = " << dt << std::endl
		  << "externFieldx = " << externFieldx << std::endl
		  << "externFieldy = " << externFieldy << std::endl
		  << "externFieldz = " << externFieldz << std::endl;
}

void Script::print(std::ofstream &ofs) {
	ofs << "# -simulation " << writeInterval << " " << timesteps << " " << dt << std::endl
		<< "#	        writeInterval timesteps dt" << std::endl
		<< "# -rectang " << nx << " " << ny << " " << nz << std::endl
        << "#        nx  ny  nz" << std::endl
		<< "# -material " << alpha << " " << exchConstant << " " << initMx << " " << initMy << " " << initMz 
        << " " << Hkx << " " << Hky << " " << Hkz << std::endl
        << "#         alpha A (J/m) M_init.x M_init.y M_init.z Hkx Hky Hkz" << std::endl
		<< "# -globalfield " << externFieldx << " " << externFieldy << " " << externFieldz << std::endl
        << "# globalField.x  y  z" << std::endl;
	/*ofs   << "# " << "nx = " << nx << std::endl
          << "# " << "ny = " << ny << std::endl
          << "# " << "nz = " << nz << std::endl
		  << "# " << "Ms = " << Ms << std::endl
		  << "# " << "alpha = " << alpha << std::endl
		  << "# " << "exchConstant = " << exchConstant << std::endl
		  << "# " << "initMx = " << initMx << std::endl
		  << "# " << "initMy = " << initMy << std::endl
		  << "# " << "initMz = " << initMz << std::endl
		  << "# " << "Hkx = " << Hkx << std::endl
		  << "# " << "Hky = " << Hky << std::endl
		  << "# " << "Hkz = " << Hkz << std::endl
		  << "# " << "writeInterval = " << writeInterval << std::endl
		  << "# " << "timesteps = " << timesteps << std::endl
		  << "# " << "dt = " << dt << std::endl
		  << "# " << "externFieldx = " << externFieldx << std::endl
		  << "# " << "externFieldy = " << externFieldy << std::endl
		  << "# " << "externFieldz = " << externFieldz << std::endl;*/
}

void Script::readRunID (std::fstream& fs) {
    fs >> runID;
    std::cout << "runID = " << runID << std::endl;
    char* junk = new char[1024];
    fs.getline(junk, 1024);
}

void Script::readGlobalField(std::fstream& fs) {
	char* junk = new char[1024];
	fs >> externFieldx >> externFieldy >> externFieldz;
	fs.getline(junk, 1024);
}

void Script::readExternField(std::fstream& fs) {
	char* junk = new char[1024];
	fs >> myExtH[myExtH->numOfFields].externHx
		>> myExtH[myExtH->numOfFields].externHy
		>> myExtH[myExtH->numOfFields].externHz
		>> myExtH[myExtH->numOfFields].startTime
		>> myExtH[myExtH->numOfFields].decayTime
		>> myExtH[myExtH->numOfFields].endTime;

	myExtH->numOfFields++; //increment the static field counter
	//myExtH[myExtH->numOfFields].print();
	fs.getline(junk, 1024);
}

void Script::readRectang(std::fstream& fs) {
	char* junk = new char[1024];
	int x, y, z;
	x=y=z=0;
	fs >> x >> y >> z;
	nx = x * 2;
	ny = y * 2;
	nz = z * 2;
	fs.getline(junk, 1024);
}

void Script::readMaterial(std::fstream& fs) {
	char* junk = new char[1024];
	fs >> alpha >> exchConstant 
	   >> initMx >> initMy >> initMz
	   >> Hkx >> Hky >> Hkz;
	fs.getline(junk, 1024);
}

void Script::readSimulation(std::fstream& fs) {
	char* junk = new char[1024];
	fs >> writeInterval >> timesteps >> dt;
	fs.getline(junk, 1024);
}

void Script::readOutput(std::fstream& fs) {}