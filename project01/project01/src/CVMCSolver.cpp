#include "CVMCSolver.h"
#include "lib.h"
#include <armadillo>
#include <iostream>

using namespace arma;
using namespace std;

VMCSolver::VMCSolver() :
    nDimensions(3),
    charge(2),
    stepLength(1.0),
    nParticles(2),
    h(0.001),
    h2(1000000),
    idum(-1),
    alpha(0.5*charge),
    beta(0.5*charge),
    nCycles(1000000),
    nAccepted(0)
{
}

double VMCSolver::runMonteCarloIntegration(
        const int &nCycles_,
        const double &stepLength_,
        const double &alpha_,
        const double &beta_)
{
    rOld.zeros(nParticles, nDimensions);
    rNew.zeros(nParticles, nDimensions);

    nCycles = nCycles_;
    alpha = alpha_;
    beta = beta_;
    stepLength = stepLength_;


    double waveFunctionOld = 0, waveFunctionNew = 0;
    double energySum = 0, energySquaredSum = 0;
    double exactEnergySum = 0, exactEnergySquaredSum = 0;
    double deltaE;

    // initial trial positions
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rOld(i,j) = stepLength * (ran2(&idum) - 0.5);
        }
    }
    rNew = rOld;

    // loop over Monte Carlo cycles
    for (int cycle = 0; cycle < nCycles; cycle++) {

        // Store the current value of the wave function
        waveFunctionOld = waveFunction(rOld);

        for (int i = 0; i < nParticles; i++) {
            // New position to test
            for (int j = 0; j < nDimensions; j++) {
                rNew(i,j) = rOld(i,j) + stepLength*(ran2(&idum) - 0.5);
            }

            // Recalculate the value of the wave function
            waveFunctionNew = waveFunction(rNew);

            // Check for step acceptance (if yes, update position, if no, reset position)
            if (ran2(&idum) <= (waveFunctionNew*waveFunctionNew) / (waveFunctionOld*waveFunctionOld)) {
                rOld.row(i) = rNew.row(i); // update position
                waveFunctionOld = waveFunctionNew;
                nAccepted++;

            } else {
                rNew = rOld; // reset position, throw away the test-position
            }
            // update energies
            deltaE = localEnergy(rNew);
            energySum += deltaE;
//            energySquaredSum += deltaE*deltaE;

//            deltaE = exactLocalEnergy(rNew);
//            exactEnergySum += deltaE;
//            exactEnergySquaredSum += deltaE*deltaE;
        }
    }
    double energy = energySum/(nCycles*nParticles);
//    double energySquared = energySquaredSum/(nCycles*nParticles);
//    double exactEnergy = exactEnergySum/(nCycles*nParticles);
//    double exactEnergySquared = exactEnergySquaredSum/(nCycles*nParticles);
//    cout << "Energy: " << energy << " Energy (squared sum): " << energySquared << endl;
//    cout << "Exact energy: " << exactEnergy << " Exact energy (squared sum): " << exactEnergySquared << endl;

    cout << "Energy = " << energy << endl;
    return energy;
}

void VMCSolver::setParameters(
        const double &newStepLength,
        const double &newAlpha,
        const double &newBeta
)
{
    stepLength = newStepLength;
    alpha = newAlpha*charge;
    beta = newBeta*charge;
}

double VMCSolver::localEnergy(const mat &r)
{
    // Kinetic energy
    mat rPlus(r), rMinus(r);
    double waveFunctionMinus, waveFunctionPlus;
    double waveFunctionCurrent = waveFunction(r);
    double kineticEnergy = 0;
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rPlus(i,j) += h;
            rMinus(i,j) -= h;
            waveFunctionMinus = waveFunction(rMinus);
            waveFunctionPlus = waveFunction(rPlus);
            //kineticEnergy -= (waveFunctionMinus + waveFunctionPlus - 2 * waveFunctionCurrent);
            kineticEnergy -= waveFunctionMinus + waveFunctionPlus;
            rPlus(i,j) = rMinus(i,j) = r(i,j);
        }
    }
    kineticEnergy = 0.5*h2*(kineticEnergy/waveFunctionCurrent
                            + 2.0*double(nParticles*nDimensions));
    //kineticEnergy = 0.5 * h2 * kineticEnergy / waveFunctionCurrent;

    // Potential energy
    double potentialEnergy = 0;
    for (int i = 0; i < nParticles; i++) {
        // r(i) = norm(r.row(i), 2);
        potentialEnergy += 1.0/norm(r.row(i), 2);
    }
    potentialEnergy *= -charge;

    // Contribution from electron-electron potential
    double r12 = 0;
    for (int i = 0; i < nParticles; i++) {
        for (int j = i + 1; j < nParticles; j++) {
            r12 = 0;
            for (int k = 0; k < nDimensions; k++) {
                r12 += (r(i,k) - r(j,k))*(r(i,k) - r(j,k));
            }
            potentialEnergy += 1/sqrt(r12);
        }
    }

    return kineticEnergy + potentialEnergy;
}

double VMCSolver::exactLocalEnergy(const mat &r)
{
    double EL1, EL2, rTemp;
    double rInverseSum = 0, rSum = 0, rijSum = 0, rijInverseSum = 0;

    for (int i = 0; i < nParticles; i++) {
        rTemp = norm(r.row(i), 2);
        rSum += rTemp;
        rInverseSum += 1.0/rTemp;
        for (int j = i + 1; j < nParticles; j++) {
            rTemp = norm(r.row(i) - r.row(j), 2);
            rijSum += rTemp;
            rijInverseSum += 1.0/rTemp;
        }
    }

//    vec3 r1vec = static_cast<vec3>(r.row(0).t());
//    vec3 r2vec = static_cast<vec3>(r.row(1).t());
    vec3 r1vec(r.row(0).t());
    vec3 r2vec(r.row(1).t());

    double r1 = norm(r1vec, 2);
    double r2 = norm(r2vec, 2);

    double betafactor = 1 + beta*rijSum;
    double rfactor = 1.0 - dot(r1vec, r2vec)/r1/r2;

    EL1 = (alpha - charge)*rInverseSum + 1.0/rijSum - alpha*alpha;
    EL2 = EL1 + (1.0/(2.0*betafactor*betafactor))*
            (
                alpha*rSum*rfactor/rijSum - 1.0/2.0/(betafactor*betafactor) - 2.0/rijSum
                + 2*beta/betafactor
            );

    return EL2;
}

double VMCSolver::dr(const mat &r)
{
    return norm(rNew.row(0) - rNew.row(1), 2);
}

double VMCSolver::waveFunction2(const mat &r)
{
    double argument = 0;
    for (int i = 0; i < nParticles; i++) {
        double rSingleParticle = 0;
        for (int j = 0; j < nDimensions; j++) {
            rSingleParticle += r(i,j) * r(i,j);
        }
        argument += sqrt(rSingleParticle);
    }

    return exp(-argument * alpha);
}

double VMCSolver::waveFunction(const mat &R)
{
    double r = 0.0;
    double rSum = 0.0;

    // r1 + r2 + ...
    for (int i = 0; i < nParticles; i++) {
        r = 0.0;
        for (int j = 0; j < nDimensions; j++)
        {
            r += R(i,j)*R(i,j);
        }
        rSum += sqrt(r);
    }

    double r12Sum = 0.0;
    // r1*r2*...
    for (int i = 0; i < nParticles; i++) {
        for (int j = i + 1; j < nParticles; j++) {
            r = 0.0;
            for (int k = 0; k < nDimensions; k++)
            {
                r += (R(i,k) - R(j,k))*(R(i,k) - R(j,k));
            }
            r12Sum += sqrt(r);
        }
    }

    return exp(-alpha*rSum + r12Sum/(2.0*(1.0 + beta*r12Sum)));
}
