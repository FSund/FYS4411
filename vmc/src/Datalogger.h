#ifndef DATALOGGER_H
#define DATALOGGER_H

#include <armadillo>
#include <fstream>

using namespace std;
using namespace arma;

class Datalogger
{
public:
    Datalogger(
            const int &myRank,
            const int &numprocs,
            const bool &binary = false);
    ~Datalogger();
    void initialize(const int &nCycles, const int &nParticles, const string &filename);
    void log(const double &dataToLog);
    void writeToFile();
    void finish();
private:
    int myRank, numprocs;
    int N;
    int cycle;
    int maxCycles;
    double* data;
    ofstream ofile;
    bool binary;
};

#endif // DATALOGGER_H
