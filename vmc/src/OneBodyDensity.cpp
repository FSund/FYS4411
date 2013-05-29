#include <src/OneBodyDensity.h>

OneBodyLogger::OneBodyLogger(
        const int &myRank,
        const int &numprocs,
        const int &nParticles) :
    myRank(myRank),
    numprocs(numprocs),
    maxCycles(1e6),
    binary(true),
    nParticles(nParticles),
    nFiles(4)
{
    if (binary)
        outFiles = new ofstream[nFiles];
    else
        outFiles = new ofstream[1];
}

OneBodyLogger::~OneBodyLogger()
{
    delete[] outFiles;
}

void OneBodyLogger::initialize(const int &nCycles)
{
    N = nCycles;
    data = mat(maxCycles*nParticles, nFiles); // each process creates individual matrix
    cycle = 0;
    string fileNameBase = "onebodylog";
    ostringstream fileName;

    char ending[4];
    ending[0] = 'x';
    ending[1] = 'y';
    ending[2] = 'z';
    ending[3] = 'r';

    if (binary)
    {
        for (int i = 0; i < nFiles; i++)
        {
            fileName << fileNameBase << "_" << myRank << "of" << numprocs << "_nParticles" << nParticles << "_nSamples" << N << "_" << ending[i] << ".bin";
            outFiles[i].open(fileName.str().c_str(), ios_base::out | ios_base::binary);
            fileName.str("");
        }
    }
    else
    {
        fileName << fileNameBase << "_" << myRank << "of" << numprocs << "_nParticles" << nParticles << "_nSamples" << N << ".dat";
        outFiles[0].open(fileName.str().c_str(), ios_base::out);
    }
}

void OneBodyLogger::log(const mat &dataToLog)
{
    int k = 0;
    for (int i = cycle*nParticles; i < (cycle+1)*nParticles; i++)
    {
        r = 0.0;
        for (int j = 0; j < 3; j++)
        {
            r += dataToLog(k, j)*dataToLog(k, j);
            data(i,j) = dataToLog(k, j);
        }
        data(i,3) = sqrt(r);
        k++;
    }

    cycle++;

    if (cycle >= maxCycles)
    {
        writeToFile();
        cycle = 0;
    }
}

void OneBodyLogger::writeToFile()
{
    if (binary)
    {
//        data.save(ofile, raw_binary);
        vec datacol(maxCycles,1);
        for (int i = 0; i < nFiles; i++)
        {
            datacol = data.col(i);
            datacol.save(outFiles[i], raw_binary);
        }
    }
    else
        data.save(outFiles[0], raw_ascii);

    cout << "OneBodyLogger writing to file" << endl;
}

void OneBodyLogger::finish()
{
    data.resize(cycle*nParticles, 4);

//    cout << "hei" << endl;

//    cout << data << endl;

    if (binary)
    {
        vec datacol(cycle*nParticles,1);
        for (int i = 0; i < nFiles; i++)
        {
            datacol = data.col(i);
            datacol.save(outFiles[i], raw_binary);
        }

    }
    else
        data.save(outFiles[0], raw_ascii);

    cout << "OneBodyLogger finishing up" << endl;
}




