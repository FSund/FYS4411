#include <unittest++/UnitTest++.h>
#include <armadillo>
#include <src/Jastrow.h>
#include <src/Slater.h>
#include <src/Orbitals.h>

using namespace arma;

double findRealSlaterRatio(int nParticles, double alpha, mat rOld, mat rNew, int currentParticle);

TEST(SlaterRatioTest) {
    // test for Beryllium
    int nParticles = 4;
    int nDimensions = 3;
    double alpha = (conv_to<double>::from(randu(1,1)))*nParticles;
    int currentParticle = 0;

    mat rOld(nParticles, nDimensions);
    mat rNew(nParticles, nDimensions);
    rOld << 1 << 1 << 1 << endr
         << 2 << 3 << 4 << endr
         << 3 << 4 << 5 << endr
         << 4 << 5 << 6 << endr;
    rNew << 2 << 2 << 2 << endr
         << 2 << 3 << 4 << endr
         << 3 << 4 << 5 << endr
         << 4 << 5 << 6 << endr;

    Slater slater(nParticles);
    slater.setAlpha(alpha);
    slater.initalize(rOld);
    slater.updatePositionAndCurrentParticle(rNew, currentParticle);

    double realRatio = findRealSlaterRatio(nParticles, alpha, rOld, rNew, currentParticle);

    cout << "Slater tests" << endl;
    printf("%.20f\n", slater.getRatio());
    printf("%.20f\n", realRatio);

    CHECK(slater.getRatio() == realRatio);

//    slater.acceptMove();
    // now the ratio should be 1 again, since we haven't given Slater a new
    // position
//    CHECK(slater.getRatio() == 1);

//    mat slaterOld = zeros<mat>(nParticles, nParticles);
//    mat slaterINVold = zeros<mat>(nParticles, nParticles);
//    mat slaterNew = zeros<mat>(nParticles, nParticles);
//    mat slaterINVnew = zeros<mat>(nParticles, nParticles);
//    for (int i = 0; i < nParticles; i++) // loop over particles
//    {
//        for (int j = 0; j < nParticles/2; j++) // loop over orbitals
//        {
//            slaterOld(i,j) = orbitals.wavefunction(rOld.row(i), 0);
//            slaterOld(i,j+nParticles/2) = orbitals.wavefunction(rOld.row(i), 1);
//            slaterNew(i,j) = orbitals.wavefunction(rNew.row(i), 0);
//            slaterNew(i,j+nParticles/2) = orbitals.wavefunction(rNew.row(i), 1);
//        }
//    }
//    slaterINVold = inv(slaterOld);
//    slaterINVnew = inv(slaterNew);
//    cout << "det(old) = " << det(slaterOld) << endl;
//    cout << "det(new) = " << det(slaterNew) << endl;
//    cout << "inv(old) = " << slaterINVold;
//    cout << "inv(new) = " << slaterINVnew;
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
//    slater.initalize(rOld);


//    slater.updatePositionAndCurrentParticle(rNew, currentParticle);
//    double realRatio = findRealSlaterRatio(nParticles, alpha, rOld, rNew, currentParticle);
//    CHECK(slater.getRatio() == realRatio);
//}

TEST(JastrowRatioTest) {
    int nParticles = 4;
    int nDimensions = 3;
    double beta = (conv_to<double>::from(randu(1,1)))*nParticles;
    int currentParticle = 1;
    double realRatio = 0.0;

    mat rOld(nParticles, nDimensions);
    mat rNew(nParticles, nDimensions);
    rOld << 1 << 1 << 1 << endr
         << 2 << 3 << 4 << endr
         << 3 << 4 << 5 << endr
         << 4 << 5 << 6 << endr;
    rNew << 1 << 1 << 1 << endr
         << 2.5 << 3.5 << 4.5 << endr
         << 3 << 4 << 5 << endr
         << 4 << 5 << 6 << endr;

    Jastrow jastrow(nParticles);
    jastrow.setBeta(beta);
    jastrow.initalize(rOld);
    jastrow.updatePositionAndCurrentParticle(rNew, currentParticle);

    realRatio = jastrow.wavefunction(rNew)/jastrow.wavefunction(rOld);

    cout << "Jastrow test" << endl;
    printf("%.20f\n", jastrow.getRatio());
    printf("%.20f\n", realRatio);

    CHECK((jastrow.getRatio() - realRatio) < 1e-15);
}

int main()
{
    return UnitTest::RunAllTests();
}

double findRealSlaterRatio(int nParticles, double alpha, mat rOld, mat rNew, int currentParticle)
{
    // finding the "real" determinant and inverse determinant UP and DOWN matrices
    // using armadillo's inverse to compare with my "efficient" algorithm
    Orbitals orbitals;
    orbitals.setAlpha(alpha);
    int N = nParticles/2;
    mat slaterUPold = zeros<mat>(N, N);
    mat slaterDOWNold = zeros<mat>(N, N);
    mat slaterUPinvOld = zeros<mat>(N, N);
    mat slaterDOWNinvOld = zeros<mat>(N, N);
    mat slaterUPnew = zeros<mat>(N, N);
    mat slaterDOWNnew = zeros<mat>(N, N);
    mat slaterUPinvNew = zeros<mat>(N, N);
    mat slaterDOWNinvNew = zeros<mat>(N, N);
    for (int i = 0; i < N; i++) // loop over particles
    {
        for (int j = 0; j < N; j++) // loop over orbitals
        {
            slaterUPold(i,j)   = orbitals.wavefunction(rOld.row(i), j);
            slaterDOWNold(i,j) = orbitals.wavefunction(rOld.row(i+N), j);
            slaterUPnew(i,j)   = orbitals.wavefunction(rNew.row(i), j);
            slaterDOWNnew(i,j) = orbitals.wavefunction(rNew.row(i+N), j);
        }
    }
    slaterUPinvOld = inv(slaterUPold);
    slaterDOWNinvOld = inv(slaterDOWNold);
    slaterUPinvNew = inv(slaterUPnew);
    slaterDOWNinvNew = inv(slaterDOWNnew);

    double realRatio = 0;
    if (currentParticle < N)
    {
        int i = currentParticle;
        for (int j = 0; j < N; j++)
            realRatio += slaterUPnew(i,j)*slaterUPinvOld(j,i);
    }
    else
    {
        int i = currentParticle - N;
        for (int j = 0; j < N; j++) // loop over orbitals
            realRatio += slaterDOWNnew(i,j)*slaterDOWNinvOld(j,i);
    }

    return realRatio;
}
