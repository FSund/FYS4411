#include <src/Datalogger.h>

Datalogger::Datalogger(const int &myRank, const int &numprocs, const bool &binary):
    myRank(myRank),
    numprocs(numprocs),
    maxCycles(1e3),
    binary(binary)
{
}

Datalogger::~Datalogger()
{
    delete data;
}

void Datalogger::initialize(const int &nCycles, const int &nParticles, const string &filename)
{
    N = nCycles;
    data = new double[maxCycles]; // each process creates individual vector
    cycle = 0;

    ostringstream fileName;
    if (binary)
    {
        fileName << filename << "_" << myRank << "of" << numprocs << ".bin";
        ofile.open(fileName.str().c_str(), ios_base::out | ios_base::binary);
    }
    else
    {
        fileName << filename << "_" << myRank << "of" << numprocs << ".dat";
        ofile.open(fileName.str().c_str(), ios_base::out);
    }
}

void Datalogger::log(const double &dataToLog)
{
    data[cycle] = dataToLog;
    cycle++;

    if (cycle >= maxCycles)
    {
        writeToFile();
        cycle = 0;
    }
}

void Datalogger::writeToFile()
{
    if (binary)
    {
//        for (int i = 0; i < maxCycles; i++)
//            ofile.write(energy[i], sizeof<)
        ofile.write((char*)data, maxCycles*sizeof(double));
    }
    else
    for (int i = 0; i < maxCycles; i++)
        ofile << data[i] << endl;

}

void Datalogger::finish()
{
    if (binary)
    {
        ofile.write((char*)data, cycle*sizeof(double));
    }
    else
    for (int i = 0; i < cycle; i++)
        ofile << data[i] << endl;
}

