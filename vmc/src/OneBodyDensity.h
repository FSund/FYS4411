#ifndef ONEBODYDENSITY_H
#define ONEBODYDENSITY_H

#include <armadillo>
#include <fstream>

using namespace std;
using namespace arma;

class OneBodyLogger
{
public:
    OneBodyLogger(const int &myRank,
                  const int &numprocs,
                  const int &nParticles);
    ~OneBodyLogger();
    void initialize(const int &nCycles);
    void log(const mat &dataToLog);
    void writeToFile();
    void finish();
private:
    int myRank, numprocs;
    int maxCycles;
    bool binary;
    int nParticles;
    int nFiles;
    int N;
    int cycle;
    mat data;
    ofstream *outFiles;
    double r;
    int ix;
};

#endif // ONEBODYDENSITY_H
