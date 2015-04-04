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