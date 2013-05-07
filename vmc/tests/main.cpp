#include <unittest++/UnitTest++.h>
#include <armadillo>
#include <src/Wavefunction.h>
#include <src/Solver/Solver.h>
#include <src/Solver/SolverMCIS.h>
#include <src/Solver/SolverMCBF.h>

using namespace arma;

double findRealSlaterRatio(int nParticles, double alpha, mat rOld, mat rNew, int currentParticle);

TEST(SlaterRatioTest) {
    // test for Beryllium
    int nParticles = 4;
    int nDimensions = 3;
    double alpha = 3.8;

    int currentParticle;
    double realRatio;
    mat rOld(nParticles, nDimensions);
    mat rNew(nParticles, nDimensions);
    Slater slater(nParticles);

    rOld << 1 << 1 << 1 << endr
         << 2 << 3 << 4 << endr
         << 3 << 4 << 5 << endr
         << 4 << 5 << 6 << endr;
    slater.setAlpha(alpha);
    slater.initialize(rOld);

    currentParticle = 0;
    rNew << 1.5 << 1.5 << 1.5 << endr
         << 2 << 3 << 4 << endr
         << 3 << 4 << 5 << endr
         << 4 << 5 << 6 << endr;
    slater.updatePositionAndCurrentParticle(rNew, currentParticle);
    realRatio = slater.wavefunction(rNew)/slater.wavefunction(rOld);
    realRatio *= realRatio;

//    cout << "Slater test 1" << endl;
//    printf("%.20f\n", slater.getRatio());
//    printf("%.20f\n", realRatio);
    CHECK(abs(slater.getRatio() - realRatio) < 1e-15);

    slater.acceptMove();
    rOld = rNew;
    currentParticle = 1;
    rNew << 1.5 << 1.5 << 1.5 << endr
         << 2.5 << 3.5 << 4.5 << endr
         << 3 << 4 << 5 << endr
         << 4 << 5 << 6 << endr;
    slater.updatePositionAndCurrentParticle(rNew, currentParticle);
    realRatio = slater.wavefunction(rNew)/slater.wavefunction(rOld);
    realRatio *= realRatio;

//    cout << "Slater test 2" << endl;
//    printf("%.20f\n", slater.getRatio());
//    printf("%.20f\n", realRatio);
    CHECK(abs(slater.getRatio() - realRatio) < 1e-15);

    slater.rejectMove();
    rNew = rOld;
    currentParticle = 2;
    rNew << 1.5 << 1.5 << 1.5 << endr
         << 2 << 3 << 4 << endr
         << 3.5 << 4.5 << 5.5 << endr
         << 4 << 5 << 6 << endr;
    slater.updatePositionAndCurrentParticle(rNew, currentParticle);
//    realRatio = findRealSlaterRatio(nParticles, alpha, rOld, rNew, currentParticle);
    realRatio = slater.wavefunction(rNew)/slater.wavefunction(rOld);
    realRatio *= realRatio;

//    cout << "Slater test 3" << endl;
//    printf("%.20f\n", slater.getRatio());
//    printf("%.20f\n", realRatio);
    CHECK(abs(slater.getRatio() - realRatio) < 1e-15);
}

TEST(JastrowRatioTest) {
    int nParticles = 4;
    int nDimensions = 3;
    double beta = 0.1;
    int currentParticle;
    double realRatio = 0.0;

    mat rOld(nParticles, nDimensions);
    mat rNew(nParticles, nDimensions);
    rOld << 1 << 1 << 1 << endr
         << 2 << 3 << 4 << endr
         << 3 << 4 << 5 << endr
         << 4 << 5 << 6 << endr;
    Jastrow jastrow(nParticles);
    jastrow.setBeta(beta);
    jastrow.initialize(rOld);

    currentParticle = 1;
    rNew << 1 << 1 << 1 << endr
         << 2.5 << 3.5 << 4.5 << endr
         << 3 << 4 << 5 << endr
         << 4 << 5 << 6 << endr;
    jastrow.updatePositionAndCurrentParticle(rNew, currentParticle);
    realRatio = jastrow.wavefunction(rNew)/jastrow.wavefunction(rOld);
    realRatio *= realRatio;
//    cout << "Jastrow test 1" << endl;
//    printf("%.20f\n", jastrow.getRatio());
//    printf("%.20f\n", realRatio);
    CHECK(abs(jastrow.getRatio() - realRatio) < 1e-15);

    jastrow.acceptMove();
    currentParticle = 2;
    rOld = rNew;
    rNew << 1 << 1 << 1 << endr
         << 2.5 << 3.5 << 4.5 << endr
         << 3.5 << 4.5 << 5.5 << endr
         << 4 << 5 << 6 << endr;
    jastrow.updatePositionAndCurrentParticle(rNew, currentParticle);
    realRatio = jastrow.wavefunction(rNew)/jastrow.wavefunction(rOld);
    realRatio *= realRatio;

//    cout << "Jastrow test 2" << endl;
//    cout << jastrow.getRatio() << endl;
//    cout << realRatio << endl;
    CHECK(abs(jastrow.getRatio() - realRatio) < 1e-15);

    jastrow.rejectMove();
    currentParticle = 3;
    rNew << 1 << 1 << 1 << endr
         << 2.5 << 3.5 << 4.5 << endr
         << 3 << 4 << 5 << endr
         << 4.5 << 5.5 << 6.5 << endr;
    jastrow.updatePositionAndCurrentParticle(rNew, currentParticle);
    realRatio = jastrow.wavefunction(rNew)/jastrow.wavefunction(rOld);
    realRatio *= realRatio;

//    cout << "Jastrow test 3" << endl;
//    cout << jastrow.getRatio() << endl;
//    cout << realRatio << endl;
    CHECK(abs(jastrow.getRatio() - realRatio) < 1e-15);
}

TEST(GradientTestNeon)
{
    int nParticles = 10;
    int nDimensions = 3;
    double alpha = 10.6;
    double beta = 0.1;
//    double h = 1e-3;
//    double tolerance = 1e-4;
    double relativeTolerance = 1e-4;

    mat gradientClosedform(nParticles, nDimensions);
    mat gradientNumerical(nParticles, nDimensions);
    mat difference(nParticles, nDimensions);
    mat rOld(nParticles, nDimensions);
    mat rNew(nParticles, nDimensions);
    mat reldiff(difference);

    rNew = rOld = randu(nParticles, nDimensions);
    Wavefunction wf(nParticles, nParticles);
    wf.setAlpha(alpha);
    wf.setBeta(beta);
    wf.initialize(rOld);

    for (int ii = 0; ii < 4; ii++)
    {
        for (int cp = 0; cp < nParticles; cp++)
        {
            rNew.row(cp) = rNew.row(cp) + randn<rowvec>(1,nDimensions);
            wf.updatePositionAndCurrentParticle(rNew, cp);

            gradientClosedform = wf.localGradient();
            gradientNumerical = wf.localGradientNumerical();

//            difference = gradientClosedform - gradientNumerical;
            reldiff = abs(gradientClosedform/gradientNumerical - 1.0);
//            cout << "reldiff" << endl << reldiff;
//            cout << "difference" << endl << difference;
            CHECK(reldiff.max() < relativeTolerance);

            if (conv_to<double>::from(randn(1,1)) < 0.5)
            {
                wf.acceptMove();
                rOld.row(cp) = rNew.row(cp);
            }
            else
            {
                wf.rejectMove();
                rNew.row(cp) = rOld.row(cp);
            }
        }
    }
}

TEST(LaplacianNeonTest)
{
    int nParticles = 10;
    int nDimensions = 3;
    double alpha = 10.6;
    double beta = 0.1;
//    double h = 1e-3;
//    double tolerance = 1e-2;
    double relativeTolerance = 1e-5;
//    double difference;
    double reldiff;
    mat rOld(nParticles, nDimensions);
    mat rNew(nParticles, nDimensions);
    double numericalLaplacian;
    double closedFormLaplacian;
    Wavefunction wf(nParticles, nParticles);
    wf.setAlpha(alpha);
    wf.setBeta(beta);
    rNew = rOld = randu(nParticles, nDimensions);
    wf.initialize(rOld);

//    cout << "rNew" << endl << rNew;
    for (int ii = 0; ii < 4; ii++)
    {
        for (int cp = 0; cp < nParticles; cp++)
        {
            rNew.row(cp) = rNew.row(cp) + randn<rowvec>(1,3);
            wf.updatePositionAndCurrentParticle(rNew, cp);

            numericalLaplacian = wf.localLaplacianNumerical();
            closedFormLaplacian = wf.localLaplacian();

//            difference = abs(numericalLaplacian - closedFormLaplacian);
            reldiff = abs(closedFormLaplacian/numericalLaplacian - 1.0);
//            cout << "difference = " << difference << endl;
//            cout << "reldiff = " << reldiff << endl;
//            CHECK(difference < tolerance);
            CHECK(reldiff < relativeTolerance);

            if (conv_to<double>::from(randn(1,1)) < 0.5)
            {
                wf.acceptMove();
                rOld.row(cp) = rNew.row(cp);
            }
            else
            {
                wf.rejectMove();
                rNew.row(cp) = rOld.row(cp);
            }
        }
    }
}

TEST(SlaterGradientNeonTest)
{
    int nParticles = 10;
    int nDimensions = 3;
    double alpha = 10.6;
    double h = 1e-3;
//    double tolerance = 1e-4;
    double relativeTolerance = 1e-4;

    mat gradientClosedform(nParticles, nDimensions);
    mat gradientNumerical(nParticles, nDimensions);
    mat difference(nParticles, nDimensions);
    mat rOld(nParticles, nDimensions);
    mat rNew(nParticles, nDimensions);
    mat reldiff(difference);

    rNew = rOld = randu(nParticles, nDimensions);
    Slater wf(nParticles);
    wf.setAlpha(alpha);
    wf.initialize(rOld);

    for (int ii = 0; ii < 4; ii++)
    {
        for (int cp = 0; cp < nParticles; cp++)
        {
            rNew.row(cp) = rNew.row(cp) + randn<rowvec>(1,nDimensions);
            wf.updatePositionAndCurrentParticle(rNew, cp);

            for (int i = 0; i < nParticles; i++)
            {
                gradientClosedform.row(i) = wf.localGradient(i);
                gradientNumerical.row(i) = wf.localGradientNumerical(i, h);
            }
            difference = gradientClosedform - gradientNumerical;
            reldiff = abs(gradientClosedform/gradientNumerical - 1.0);
//            cout << "reldiff" << endl << reldiff;
//            cout << "difference" << endl << difference;
            CHECK(reldiff.max() < relativeTolerance);

            if (conv_to<double>::from(randn(1,1)) < 0.5)
            {
                wf.acceptMove();
                rOld.row(cp) = rNew.row(cp);
            }
            else
            {
                wf.rejectMove();
                rNew.row(cp) = rOld.row(cp);
            }
        }
    }
}

TEST(JastrowGradientNeonTest)
{
    int nParticles = 10;
    int nDimensions = 3;
    double beta = 0.1;
    double h = 1e-3;
//    double tolerance = 1e-4;
    double relativeTolerance = 1e-4;

    mat gradientClosedform(nParticles, nDimensions);
    mat gradientNumerical(nParticles, nDimensions);
    mat difference(nParticles, nDimensions);
    mat reldiff(difference);
    mat rOld(nParticles, nDimensions);
    mat rNew(nParticles, nDimensions);
    Jastrow wf(nParticles);
    wf.setBeta(beta);

    rNew = rOld = randu(nParticles, nDimensions);
    wf.initialize(rOld);

    for (int ii = 0; ii < 4; ii++)
    {
        for (int cp = 0; cp < nParticles; cp++)
        {
            rNew.row(cp) = rNew.row(cp) + randn<rowvec>(1,3);
            wf.updatePositionAndCurrentParticle(rNew, cp);

            for (int i = 0; i < nParticles; i++)
            {
                gradientClosedform.row(i) = wf.localGradient(i);
                gradientNumerical.row(i) = wf.localGradientNumerical(i, h);
            }
//            difference = abs(gradientClosedform - gradientNumerical);
            reldiff = abs(gradientClosedform/gradientNumerical - 1.0);
            // cout << "difference = " << difference << endl;
//            CHECK(difference.max() < tolerance);
            CHECK(reldiff.max() < relativeTolerance);

            if (conv_to<double>::from(randn(1,1)) < 0.5)
            {
                wf.acceptMove();
                rOld.row(cp) = rNew.row(cp);
            }
            else
            {
                wf.rejectMove();
                rNew.row(cp) = rOld.row(cp);
            }
        }
    }
}

TEST(SlaterLaplacianNeonTest)
{
    int nParticles = 10;
    int nDimensions = 3;
    double alpha = 10.6;
    double h = 1e-3;
//    double tolerance = 1e-2;
    double relativeTolerance = 1e-4;
//    double difference;
    double reldiff;
    mat rOld(nParticles, nDimensions);
    mat rNew(nParticles, nDimensions);
    double numericalLaplacian;
    double closedFormLaplacian;
    Slater wf(nParticles);
    wf.setAlpha(alpha);
    rNew = rOld = randu(nParticles, nDimensions);
    wf.initialize(rOld);

//    cout << "rNew" << endl << rNew;
    for (int ii = 0; ii < 4; ii++)
    {
        for (int cp = 0; cp < nParticles; cp++)
        {
            rNew.row(cp) = rNew.row(cp) + randn<rowvec>(1,3);
            wf.updatePositionAndCurrentParticle(rNew, cp);

            numericalLaplacian = closedFormLaplacian = 0.0;
            for (int i = 0; i < nParticles; i++)
            {
                numericalLaplacian += wf.localLaplacianNumerical(i, h);
                closedFormLaplacian += wf.localLaplacian(i);
            }
//            difference = abs(numericalLaplacian - closedFormLaplacian);
            reldiff = abs(numericalLaplacian/closedFormLaplacian - 1.0);
//            cout << "difference = " << difference << endl;
//            CHECK(difference < tolerance);
            CHECK(reldiff < relativeTolerance);

            if (conv_to<double>::from(randn(1,1)) < 0.5)
            {
                wf.acceptMove();
                rOld.row(cp) = rNew.row(cp);
            }
            else
            {
                wf.rejectMove();
                rNew.row(cp) = rOld.row(cp);
            }
        }
    }
}

TEST(JastrowLaplacianNeonTest)
{
    int nParticles = 10;
    int nDimensions = 3;
//    double alpha = 10.6;
    double beta = 0.1;
    double h = 1e-3;
//    double tolerance = 1e-3;
    double relativeTolerance = 1e-5;

//    double difference;
    double reldiff;
    mat rOld(nParticles, nDimensions);
    mat rNew(nParticles, nDimensions);
    double numericalLaplacian;
    double closedFormLaplacian;

    rNew = rOld = randu(nParticles, nDimensions);
    Jastrow wf(nParticles);
    wf.setBeta(beta);
    wf.initialize(rOld);

    for (int ii = 0; ii < 4; ii++)
    {
        for (int cp = 0; cp < nParticles; cp++)
        {
            rNew.row(cp) = rNew.row(cp) + randn<rowvec>(1,3);
            wf.updatePositionAndCurrentParticle(rNew, cp);

            numericalLaplacian = closedFormLaplacian = 0.0;
            for (int i = 0; i < nParticles; i++)
            {
                numericalLaplacian += wf.localLaplacianNumerical(i, h);
                closedFormLaplacian += wf.localLaplacian(i);
            }
//            difference = abs(numericalLaplacian - closedFormLaplacian);
            reldiff = abs(numericalLaplacian/closedFormLaplacian - 1.0);
//            cout << "difference = " << difference << endl;
//            CHECK(difference < tolerance);
            CHECK(reldiff < relativeTolerance);

            if (conv_to<double>::from(randn(1,1)) < 0.5)
            {
                wf.acceptMove();
                rOld.row(cp) = rNew.row(cp);
            }
            else
            {
                wf.rejectMove();
                rNew.row(cp) = rOld.row(cp);
            }
        }
    }
}

TEST(SlaterInverseTest)
{
    int nParticles = 4;
    int nDimensions = 3;
    double alpha = 3.6;
    int currentParticle;
    mat difference;
    mat rOld(nParticles, nDimensions);
    mat rNew(nParticles, nDimensions);
    Slater wf(nParticles);

    rOld << 1 << 2 << 3 << endr
         << 2 << 3 << 4 << endr
         << 3 << 4 << 5 << endr
         << 4 << 5 << 6 << endr;
    wf.setAlpha(alpha);
    wf.initialize(rOld);

    currentParticle = 0;
    rNew << 1.5 << 2.5 << 3.5 << endr
         << 2 << 3 << 4 << endr
         << 3 << 4 << 5 << endr
         << 4 << 5 << 6 << endr;
    wf.updatePositionAndCurrentParticle(rNew, currentParticle);

    difference = abs(wf.getUPinvOld() - inv(wf.getUPold()));
    CHECK(difference.max() < 1e-4); // ok

    wf.acceptMove(); // calculating the new inverse in Slater
    rOld = rNew;
    difference = abs(wf.getUPinvOld() - inv(wf.getUPold()));
//    cout << "UPinvOld" << endl << wf.getUPinvOld();
//    cout << "inv(wf.getUPold())" << endl << inv(wf.getUPold());
//    cout << "difference.max() = " << difference.max() << endl;
    CHECK(difference.max() < 1e-4); // FAILS

    currentParticle = 2;
    rNew << 1.5 << 2.5 << 3.5 << endr
         << 2 << 3 << 4 << endr
         << 3.5 << 4.5 << 5.5 << endr
         << 4 << 5 << 6 << endr;

    wf.updatePositionAndCurrentParticle(rNew, currentParticle);

    difference = abs(wf.getUPinvOld() - inv(wf.getUPold()));
    CHECK(difference.max() < 1e-4);

    wf.rejectMove(); // calculating the new inverse in Slater
    rOld = rNew;
    difference = abs(wf.getUPinvOld() - inv(wf.getUPold()));
    CHECK(difference.max() < 1e-4);

    currentParticle = 3;
    rNew << 1.5 << 2.5 << 3.5 << endr
         << 2 << 3 << 4 << endr
         << 3 << 4 << 5 << endr
         << 4.5 << 5.5 << 6.5 << endr;

    wf.updatePositionAndCurrentParticle(rNew, currentParticle);
    wf.getRatio(); // need this to update the ratio internally in Slater

    difference = abs(wf.getUPinvOld() - inv(wf.getUPold()));
    CHECK(difference.max() < 1e-4);

    wf.acceptMove(); // calculating the new inverse in Slater
    rOld = rNew;
    difference = abs(wf.getUPinvOld() - inv(wf.getUPold()));
    CHECK(difference.max() < 1e-4);
}

//TEST(ImportanceSampling)
//{
////    double alpha = 1.8;
////    double beta = 0.36;
////    int nParticles = 2;
////    int charge = 2;

//    double alpha = 3.9;
//    double beta = 4.0;
//    int nParticles = 4;
//    int charge = 4;
//    int myRank = 0;
//    int numprocs = 1;

//    SolverMCIS solver(myRank, numprocs, nParticles, charge);
//    solver.setAlpha(alpha);
//    solver.setBeta(beta);

//    double energy = solver.testSolver(1);
//    CHECK(energy == 1);
//}

//TEST(ClosedFormEnergyTestIS)
//{
//    double alpha, beta;
//    int nParticles, charge;
//    int nCycles;
//    double energyNUM, energyCF;
//    double difference;
//    double tolerance;

//    nCycles = 1e4;
//    tolerance = 1e-4;

//    // neon
//    alpha = 10.6;
//    beta = 0.1;
//    nParticles = 10;
//    charge = nParticles;

//    SolverMCIS solver(myRank, numprocs, nParticles, charge);
//    solver.setAlpha(alpha);
//    solver.setBeta(beta);
//    solver.setClosedform(false);
//    energyNUM = solver.runMonteCarloIntegration(nCycles);
//    solver.setClosedform(true);
//    energyCF = solver.runMonteCarloIntegration(nCycles);

//    difference = abs(energyNUM - energyCF);
//    CHECK(difference < tolerance);
//}

//TEST(ClosedFormEnergyTestBF)
//{
//    double alpha, beta;
//    int nParticles, charge;
//    int nCycles;
//    double energyNUM, energyCF;
//    double difference;
//    double tolerance;

//    nCycles = 1e4;
//    tolerance = 1e-4;

//    // neon
//    alpha = 10.6;
//    beta = 0.1;
//    nParticles = 10;
//    charge = nParticles;
//    SolverMCBF solver(myRank, numprocs, nParticles, charge);
//    solver.setAlpha(alpha);
//    solver.setBeta(beta);
//    solver.setClosedform(false);
//    energyNUM = solver.runMonteCarloIntegration(nCycles);
//    solver.setClosedform(true);
//    energyCF = solver.runMonteCarloIntegration(nCycles);

//    difference = abs(energyNUM - energyCF);
//    CHECK(difference < tolerance);
//}

int main(int argc, char *argv[])
{
    // MPI initialization
    int myRank, numprocs;
    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &myRank);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);

    return UnitTest::RunAllTests();
    MPI_Finalize();
}

