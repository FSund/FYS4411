#include <unittest++/UnitTest++.h>
#include <armadillo>
#include <src/Jastrow.h>
#include <src/Slater.h>
#include <src/Orbitals.h>
#include <src/Wavefunction.h>
#include <src/Solver/Solver.h>
#include <src/Solver/SolverMCIS.h>

using namespace arma;

double findRealSlaterRatio(int nParticles, double alpha, mat rOld, mat rNew, int currentParticle);

TEST(SlaterRatioTest) {
    // test for Beryllium
    double realRatio = 0.0;
    int nParticles = 4;
    int nDimensions = 3;
    double alpha = (conv_to<double>::from(randu(1,1)))*nParticles;
    int currentParticle;

    mat rOld(nParticles, nDimensions);
    mat rNew(nParticles, nDimensions);
    rOld << 1 << 1 << 1 << endr
         << 2 << 3 << 4 << endr
         << 3 << 4 << 5 << endr
         << 4 << 5 << 6 << endr;
    Slater slater(nParticles);
    slater.setAlpha(alpha);
    slater.initialize(rOld);

    currentParticle = 0;
    rNew << 1.5 << 1.5 << 1.5 << endr
         << 2 << 3 << 4 << endr
         << 3 << 4 << 5 << endr
         << 4 << 5 << 6 << endr;
    slater.updatePositionAndCurrentParticle(rNew, currentParticle);
//    realRatio = findRealSlaterRatio(nParticles, alpha, rOld, rNew, currentParticle);
    realRatio = slater.wavefunction(rNew)/slater.wavefunction(rOld);

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
//    realRatio = findRealSlaterRatio(nParticles, alpha, rOld, rNew, currentParticle);
    realRatio = slater.wavefunction(rNew)/slater.wavefunction(rOld);

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

//    cout << "Slater test 3" << endl;
//    printf("%.20f\n", slater.getRatio());
//    printf("%.20f\n", realRatio);
    CHECK(abs(slater.getRatio() - realRatio) < 1e-15);
}

//TEST(SlaterInverseTest)
//{
//    // test for Beryllium
//    int nParticles = 4;
//    int nDimensions = 3;
//    double alpha = (conv_to<double>::from(randu(1,1)))*nParticles;
//    int currentParticle = 0;

//    mat rOld(nParticles, nDimensions);
//    mat rNew(nParticles, nDimensions);
//    rOld << 1 << 1 << 1 << endr
//         << 2 << 3 << 4 << endr
//         << 3 << 4 << 5 << endr
//         << 4 << 5 << 6 << endr;
//    rNew << 2 << 2 << 2 << endr
//         << 2 << 3 << 4 << endr
//         << 3 << 4 << 5 << endr
//         << 4 << 5 << 6 << endr;
//    rOld = randu(nParticles, nDimensions);
//    rNew = randu(nParticles, nDimensions);

//    Slater slater(nParticles);
//    slater.setAlpha(alpha);
//    slater.initialize(rOld);


//    slater.updatePositionAndCurrentParticle(rNew, currentParticle);
//    double realRatio = findRealSlaterRatio(nParticles, alpha, rOld, rNew, currentParticle);
//    CHECK(slater.getRatio() == realRatio);
//}

TEST(JastrowRatioTest) {
    int nParticles = 4;
    int nDimensions = 3;
    double beta = (conv_to<double>::from(randu(1,1)))*nParticles;
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

//    cout << "Jastrow test 3" << endl;
//    cout << jastrow.getRatio() << endl;
//    cout << realRatio << endl;
    CHECK(abs(jastrow.getRatio() - realRatio) < 1e-15);
}

TEST(GradientTest)
{
    int nParticles = 4;
    int nDimensions = 3;
    vec parameters = randu(2,1);
    int currentParticle;
    double charge = 4.0;
    mat difference;

    mat rOld(nParticles, nDimensions);
    mat rNew(nParticles, nDimensions);
    rOld << 1 << 1 << 1 << endr
         << 2 << 3 << 4 << endr
         << 3 << 4 << 5 << endr
         << 4 << 5 << 6 << endr;
    Wavefunction wf(nParticles, charge);
    wf.setParameters(parameters);
    wf.initialize(rOld);

    difference = abs(wf.localGradientNumerical() - wf.localGradient());
    CHECK(difference.max() < 1e-4);

    currentParticle = 1;
    rNew << 1 << 1 << 1 << endr
         << 2.5 << 3.5 << 4.5 << endr
         << 3 << 4 << 5 << endr
         << 4 << 5 << 6 << endr;
    wf.updatePositionAndCurrentParticle(rNew, currentParticle);
    wf.getRatio(); // to update the ratio, needed to update the inverse later
    wf.acceptMove();

    difference = abs(wf.localGradientNumerical() - wf.localGradient());
    CHECK(difference.max() < 1e-4);
}

TEST(LaplacianTest)
{
    int nParticles = 4;
    int nDimensions = 3;
    vec parameters = randu(2,1);
//    int currentParticle;
    double charge = 4.0;
    double difference;

    mat rOld(nParticles, nDimensions);
    mat rNew(nParticles, nDimensions);
    rOld << 1 << 2 << 3 << endr
         << 2 << 3 << 4 << endr
         << 3 << 4 << 5 << endr
         << 4 << 5 << 6 << endr;
    Wavefunction wf(nParticles, charge);
    wf.setParameters(parameters);
    wf.initialize(rOld);

//    cout << "lapl numerical = " << wf.localLaplacianNumerical() << endl;
//    cout << "lapl closedform = " << wf.localLaplacian() << endl;

    difference = abs(wf.localLaplacianNumerical() - wf.localLaplacian());
    CHECK(difference < 1e-4);
}

//TEST(SlaterGradientTest)
//{
//cout << "-------------------------- SlaterGradientTest --------------------------" << endl;

//    int nParticles = 4;
//    int nDimensions = 3;

//    double alpha = 3.8;
//    double h = 1e-3;

//    mat rOld(nParticles, nDimensions);
//    mat rNew(nParticles, nDimensions);
//    mat gradientClosedform(nParticles, nDimensions);
//    mat difference;
//    int currentParticle = 0;

//    Slater wf(nParticles);
//    rOld << 1 << 2 << 3 << endr
//         << 2 << 3 << 4 << endr
//         << 3 << 4 << 5 << endr
//         << 4 << 5 << 6 << endr;
//    wf.setAlpha(alpha);
//    wf.initialize(rOld);

//    ////////////////////////////////////////////
//    cout << "gradient closed form" << endl;
//    for (int i = 0; i < nParticles; i++)
//        cout << wf.localGradient(i);
//    cout << "gradient numerical" << endl;
//    for (int i = 0; i < nParticles; i++)
//        cout << wf.localGradientNumerical(i, h);
//    cout << endl;
//    ////////////////////////////////////////////

//    difference = abs(wf.localGradientNumerical(currentParticle, h) - wf.localGradient(currentParticle));
//    CHECK(difference.max() < 1e-4);

//    currentParticle = 1;
//    rNew << 1 << 2 << 3 << endr
//         << 2.2 << 3.2 << 4.2 << endr
//         << 3 << 4 << 5 << endr
//         << 4 << 5 << 6 << endr;
//    wf.updatePositionAndCurrentParticle(rNew, currentParticle);
//    wf.getRatio(); // to update the ratio, needed to update the inverse

//    ////////////////////////////////////////////
//    cout << "gradient closed form" << endl;
//    for (int i = 0; i < nParticles; i++)
//        cout << wf.localGradient(i);
//    cout << "gradient numerical" << endl;
//    for (int i = 0; i < nParticles; i++)
//        cout << wf.localGradientNumerical(i, h);
//    cout << endl;
//    ////////////////////////////////////////////

//    difference = abs(wf.localGradientNumerical(currentParticle, h) - wf.localGradient(currentParticle));
//    cout << "difference.max() = " << difference.max() << endl;
//    cout << difference;
//    CHECK(difference.max() < 1e-4); // fails

//cout << "-------------------------- SlaterGradientTest --------------------------" << endl;

//    // testing after accepting
//    wf.acceptMove();
//    difference = abs(wf.localGradientNumerical(currentParticle, h) - wf.localGradient(currentParticle));
//    CHECK(difference.max() < 1e-4);

//    currentParticle = 2;
//    rNew << 1 << 2 << 3 << endr
//         << 2.5 << 3.5 << 4.5 << endr
//         << 3.5 << 4.5 << 5.5 << endr
//         << 4 << 5 << 6 << endr;
//    wf.updatePositionAndCurrentParticle(rNew, currentParticle);

//    ////////////////////////////////////////////
//    cout << "gradient closed form" << endl;
//    for (int i = 0; i < nParticles; i++)
//        cout << wf.localGradient(i);
//    cout << "gradient numerical" << endl;
//    for (int i = 0; i < nParticles; i++)
//        cout << wf.localGradientNumerical(i, h);
//    cout << endl;
//    ////////////////////////////////////////////

//    // testing before rejecting (for quantum force)
//    difference = abs(wf.localGradientNumerical(currentParticle, h) - wf.localGradient(currentParticle));
//    cout << "difference.max() = " << difference.max() << endl;
//    cout << difference;
//    CHECK(difference.max() < 1e-4); // FAILS!

//    // testing after rejecting
//    wf.rejectMove();
//    difference = abs(wf.localGradientNumerical(currentParticle, h) - wf.localGradient(currentParticle));
//    CHECK(difference.max() < 1e-4);
//}

//TEST(SlaterGradientTest)
//{
//cout << "-------------------------- SlaterGradientTest --------------------------" << endl;

//    int nParticles = 4;
//    int nDimensions = 3;

//    double alpha = 3.8;
//    double h = 1e-3;

//    mat rOld(nParticles, nDimensions);
//    mat rNew(nParticles, nDimensions);
//    mat gradientClosedform(nParticles, nDimensions);
//    mat difference;
//    int currentParticle = 0;

//    Slater wf(nParticles);
//    rOld << 1 << 2 << 3 << endr
//         << 2 << 3 << 4 << endr
//         << 3 << 4 << 5 << endr
//         << 4 << 5 << 6 << endr;
//    wf.setAlpha(alpha);
//    wf.initialize(rOld);

//    ////////////////////////////////////////////
//    cout << "gradient closed form" << endl;
//    for (int i = 0; i < nParticles; i++)
//        cout << wf.localGradient(i);
//    cout << "gradient closed form2" << endl;
//    for (int i = 0; i < nParticles; i++)
//        cout << wf.localGradient(rOld, i);
//    cout << "gradient numerical" << endl;
//    for (int i = 0; i < nParticles; i++)
//        cout << wf.localGradientNumerical(i, h);
//    cout << "gradient numerical2" << endl;
//    cout << wf.gradient(h);
//    cout << "gradient numerical3" << endl;
//    cout << wf.gradient(rOld, h);
//    cout << endl;
//    ////////////////////////////////////////////

//    difference = abs(wf.localGradientNumerical(currentParticle, h) - wf.localGradient(currentParticle));
//    CHECK(difference.max() < 1e-4);

//    currentParticle = 1;
//    rNew << 1 << 2 << 3 << endr
//         << 2.2 << 3.2 << 4.2 << endr
//         << 3 << 4 << 5 << endr
//         << 4 << 5 << 6 << endr;
//    wf.updatePositionAndCurrentParticle(rNew, currentParticle);
//    wf.getRatio(); // to update the ratio, needed to update the inverse

//    ////////////////////////////////////////////
//    cout << "gradient closed form" << endl;
//    for (int i = 0; i < nParticles; i++)
//        cout << wf.localGradient(i);
//    cout << "gradient closed form2" << endl;
//    for (int i = 0; i < nParticles; i++)
//        cout << wf.localGradient(rNew, i);
//    cout << "gradient numerical" << endl;
//    for (int i = 0; i < nParticles; i++)
//        cout << wf.localGradientNumerical(i, h);
//    cout << "gradient numerical2" << endl;
//    cout << wf.gradient(h);
//    cout << "gradient numerical3" << endl;
//    cout << wf.gradient(rNew, h);
//    cout << endl;
//    ////////////////////////////////////////////

//    difference = abs(wf.localGradientNumerical(currentParticle, h) - wf.localGradient(currentParticle));
//    cout << "difference.max() = " << difference.max() << endl;
//    cout << difference;
//    CHECK(difference.max() < 1e-4); // fails

//cout << "-------------------------- SlaterGradientTest --------------------------" << endl;
//}

TEST(JastrowGradientTest)
{
    int nParticles = 4;
    int nDimensions = 3;

    double beta = 0.1;
    double h = 1e-3;

    int currentParticle;
    mat gradientClosedform = zeros<mat>(nParticles, nDimensions);

    mat rOld(nParticles, nDimensions);
    mat rNew(nParticles, nDimensions);
    rOld << 1 << 2 << 3 << endr
         << 2 << 3 << 4 << endr
         << 3 << 4 << 5 << endr
         << 4 << 5 << 6 << endr;
    Jastrow wf(nParticles);
    wf.setBeta(beta);
    wf.initialize(rOld);

    mat difference;
    difference = abs(wf.localGradientNumerical(currentParticle, h) - wf.localGradient(currentParticle));
    CHECK(difference.max() < 1e-4);

    currentParticle = 1;
    rNew << 1 << 2 << 3 << endr
         << 2.5 << 3.5 << 4.5 << endr
         << 3 << 4 << 5 << endr
         << 4 << 5 << 6 << endr;
    wf.updatePositionAndCurrentParticle(rNew, currentParticle);

    // test before accepting (because the quantum force uses the gradient
    // before accepting/rejecting
    difference = abs(wf.localGradientNumerical(currentParticle, h) - wf.localGradient(currentParticle));
    CHECK(difference.max() < 1e-4);

    // testing after accepting
    wf.acceptMove();
    for (int i = 0; i < nParticles; i++)
        gradientClosedform.row(i) = wf.localGradient(i);
    difference = abs(wf.localGradientNumerical(currentParticle, h) - wf.localGradient(currentParticle));
    CHECK(difference.max() < 1e-4);

    currentParticle = 2;
    rNew << 1 << 2 << 3 << endr
         << 2.5 << 3.5 << 4.5 << endr
         << 3.5 << 4.5 << 5.5 << endr
         << 4 << 5 << 6 << endr;
    wf.updatePositionAndCurrentParticle(rNew, currentParticle);

    // testing before rejecting (for quantum force)
    difference = abs(wf.localGradientNumerical(currentParticle, h) - wf.localGradient(currentParticle));
    CHECK(difference.max() < 1e-4);

    // testing after rejecting
    wf.rejectMove();
    difference = abs(wf.localGradientNumerical(currentParticle, h) - wf.localGradient(currentParticle));
    CHECK(difference.max() < 1e-4);
}

TEST(JastrowLaplacianTest)
{
    int nParticles = 4;
    int nDimensions = 3;
    double beta = 0.1;
//    int currentParticle;
//    double charge = 4.0;
    double difference;
    double h = 1e-3;
    int currentParticle = 0;

    mat rOld(nParticles, nDimensions);
    mat rNew(nParticles, nDimensions);
    rOld << 1 << 2 << 3 << endr
         << 2 << 3 << 4 << endr
         << 3 << 4 << 5 << endr
         << 4 << 5 << 6 << endr;
    Jastrow wf(nParticles);
    wf.setBeta(beta);
    wf.initialize(rOld);

//    cout << "lapl numerical = " << jastrow.localLaplacianNumerical(h) << endl;
//    cout << "lapl closedform = " << laplacianClosedform << endl;

    difference = abs(wf.localLaplacianNumerical(currentParticle, h) - wf.localLaplacian(currentParticle));
    CHECK(difference < 1e-4);
}

TEST(SlaterLaplacianTest)
{
    int nParticles = 4;
    int nDimensions = 3;
    double alpha = 3.6;
//    int currentParticle;
//    double charge = 4.0;
    double difference;
    double h = 1e-3;
    int currentParticle = 0;
    mat rOld(nParticles, nDimensions);
    mat rNew(nParticles, nDimensions);

    rOld << 1 << 2 << 3 << endr
         << 2 << 3 << 4 << endr
         << 3 << 4 << 5 << endr
         << 4 << 5 << 6 << endr;
    Slater wf(nParticles);
    wf.setAlpha(alpha);
    wf.initialize(rOld);

//    cout << "lapl numerical = " << slater.localLaplacianNumerical(h) << endl;
//    cout << "lapl closedform = " << laplacianClosedform << endl;

    difference = abs(wf.localLaplacianNumerical(currentParticle, h) - wf.localLaplacian(currentParticle));
    CHECK(difference < 1e-4);
}

TEST(SlaterInverseTest)
{
    int nParticles = 4;
    int nDimensions = 3;
    double alpha = 3.6;
    int currentParticle;
//    double charge = 4.0;
    mat difference;
//    double h = 1e-3;
    mat rOld(nParticles, nDimensions);
    mat rNew(nParticles, nDimensions);

    rOld << 1 << 2 << 3 << endr
         << 2 << 3 << 4 << endr
         << 3 << 4 << 5 << endr
         << 4 << 5 << 6 << endr;
    Slater wf(nParticles);
    wf.setAlpha(alpha);
    wf.initialize(rOld);

    currentParticle = 0;
    rNew << 1.5 << 2.5 << 3.5 << endr
         << 2 << 3 << 4 << endr
         << 3 << 4 << 5 << endr
         << 4 << 5 << 6 << endr;
    wf.updatePositionAndCurrentParticle(rNew, currentParticle);
    wf.getRatio(); // need this to update the ratio internally in slater

    difference = abs(wf.getUPinvOld() - inv(wf.getUPold()));
    CHECK(difference.max() < 1e-4);

    wf.acceptMove();

    difference = abs(wf.getUPinvNew() - inv(wf.getUPnew()));
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

//double findRealSlaterRatio(int nParticles, double alpha, mat rOld, mat rNew, int currentParticle)
//{
//    // finding the "real" determinant and inverse determinant UP and DOWN matrices
//    // using armadillo's inverse to compare with my "efficient" algorithm
//    Orbitals orbitals;
//    orbitals.setAlpha(alpha);
//    int N = nParticles/2;
//    mat slaterUPold = zeros<mat>(N, N);
//    mat slaterDOWNold = zeros<mat>(N, N);
////    mat slaterUPinvOld = zeros<mat>(N, N);
////    mat slaterDOWNinvOld = zeros<mat>(N, N);
//    mat slaterUPnew = zeros<mat>(N, N);
//    mat slaterDOWNnew = zeros<mat>(N, N);
////    mat slaterUPinvNew = zeros<mat>(N, N);
////    mat slaterDOWNinvNew = zeros<mat>(N, N);
//    for (int i = 0; i < N; i++) // loop over particles
//    {
//        for (int j = 0; j < N; j++) // loop over orbitals
//        {
//            slaterUPold(i,j)   = orbitals.wavefunction(rOld.row(i), j);
//            slaterDOWNold(i,j) = orbitals.wavefunction(rOld.row(i+N), j);
//            slaterUPnew(i,j)   = orbitals.wavefunction(rNew.row(i), j);
//            slaterDOWNnew(i,j) = orbitals.wavefunction(rNew.row(i+N), j);
//        }
//    }
////    slaterUPinvOld = inv(slaterUPold);
////    slaterDOWNinvOld = inv(slaterDOWNold);
////    slaterUPinvNew = inv(slaterUPnew);
////    slaterDOWNinvNew = inv(slaterDOWNnew);

//    double realRatio = 0;
////    if (currentParticle < N)
////    {
////        int i = currentParticle;
////        for (int j = 0; j < N; j++)
////            realRatio += slaterUPnew(i,j)*slaterUPinvOld(j,i);
////    }
////    else
////    {
////        int i = currentParticle - N;
////        for (int j = 0; j < N; j++) // loop over orbitals
////            realRatio += slaterDOWNnew(i,j)*slaterDOWNinvOld(j,i);
////    }

//    if (currentParticle < N)
//        realRatio = det(slaterUPnew)/det(slaterUPold);
//    else
//        realRatio = det(slaterDOWNnew)/det(slaterDOWNold);

//    return realRatio;
//}
