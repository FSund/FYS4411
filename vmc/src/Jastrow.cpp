#include <src/Jastrow.h>

Jastrow::Jastrow()
{
    cout << "! Error: Using default constructor in Jastrow ! " << endl;
    exit(1);
}

Jastrow::Jastrow(const int &nParticles):
    nParticles(nParticles),
    nDimensions(3),
    rOld(zeros<mat>(nParticles, nDimensions)),
    rNew(zeros<mat>(nParticles, nDimensions)),
    rijOld(zeros<mat>(nParticles, nParticles)),
    rijNew(zeros<mat>(nParticles, nParticles)),
    fijOld(zeros<mat>(nParticles, nParticles)),
    fijNew(zeros<mat>(nParticles, nParticles)),
    a(zeros<mat>(nParticles, nParticles))
{
    for (int i = 0; i < nParticles; i++)
        for (int j = 0; j < nParticles; j++)
            a(i,j) = ((i+j)%2 == 0) ? 0.25 : 0.5;
}

void Jastrow::initialize(const mat &r)
{
    rNew = rOld = r;

    calculate_rij(); // calculates rijOld from rOld
    calculate_fij(); // calculates fijOld from rOld/rijOld

    rijNew = rijOld;
    fijNew = fijOld;
}

void Jastrow::updatePositionAndCurrentParticle(mat &r, int &k)
{
    rNew = r;
    currentParticle = k;

    update_rij();
    update_fij();
}

void Jastrow::setBeta(const double &newBeta)
{
    beta = newBeta;
}

double Jastrow::getRatio()
{
    double dU = 0.0;
    int k = currentParticle;

    for (int i = 0; i < k; i++)
        dU += fijNew(i, k) - fijOld(i, k);
    for (int i = k+1; i < nParticles; i++)
        dU += fijNew(k, i) - fijOld(k, i);

    return exp(2.0*dU); // squared
}

void Jastrow::acceptMove()
{
    rOld.row(currentParticle) = rNew.row(currentParticle);

    fijOld.row(currentParticle) = fijNew.row(currentParticle);
    fijOld.col(currentParticle) = fijNew.col(currentParticle);
//    fijOld = fijNew;
    rijOld.row(currentParticle) = rijNew.row(currentParticle);
    rijOld.col(currentParticle) = rijNew.col(currentParticle);
//    rijOld = rijNew;
}

void Jastrow::rejectMove()
{
    rNew.row(currentParticle) = rOld.row(currentParticle);

    fijNew.row(currentParticle) = fijOld.row(currentParticle);
    fijNew.col(currentParticle) = fijOld.col(currentParticle);
//    fijNew = fijOld;
    rijNew.row(currentParticle) = rijOld.row(currentParticle);
    rijNew.col(currentParticle) = rijOld.col(currentParticle);
//    rijOld = rijNew;
}

double Jastrow::wavefunction()
{
    /* Only for use in the numerical gradient and laplacian in the local energy */

    double arg = 0.0;
    for (int i = 0; i < nParticles; i++)
        for (int j = i+1; j < nParticles; j++)
            arg += fijNew(i,j);

    return exp(arg);
}

double Jastrow::wavefunction(const mat &r)
{
    /* Slow function, recalculates whole rij and fij matrices!
     * Should only be used for the numerical calculation of the local energy */

    // temp rij matrix for this r-matrix
    mat rijTemp = zeros<mat>(nParticles, nParticles);
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = i+1; j < nParticles; j++)
        {
            for (int k = 0; k < nDimensions; k++)
            {
                rijTemp(i,j) += (r(i,k) - r(j,k))*(r(i,k) - r(j,k));
            }
            rijTemp(i,j) = sqrt(rijTemp(i,j));
        }
    }

    // temp fij matrix for this r-matrix
    mat fijTemp = zeros<mat>(nParticles, nParticles);
    double rij;
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < nParticles; j++)
        {
            rij = rijTemp(i,j);
            fijTemp(i,j) = a(i,j)*rij/(1.0 + beta*rij);
        }
    }

    // the wavefunction
    double arg = 0.0;
    for (int i = 0; i < nParticles; i++)
        for (int j = i+1; j < nParticles; j++)
            arg += fijTemp(i,j);

    return exp(arg);
}

rowvec Jastrow::localGradient(const int &k)
{
    /* Calculates the gradient for particle k */
    grad = zeros<rowvec>(nDimensions);
    for (int i = 0; i < k; i++)
    {
        rij = rijNew(i,k);
        grad += (rNew.row(k) - rNew.row(i)) * a(i,k)/(rij*pow((1.0 + beta*rij),2));
    }
    for (int i = k+1; i < nParticles; i++)
    {
        rij = rijNew(k,i);
        grad += (rNew.row(k) - rNew.row(i)) * a(i,k)/(rij*pow((1.0 + beta*rij),2));
    }

    return grad;
}

rowvec Jastrow::localGradientNumerical(const mat &r, const int &k, const double &h)
{
    /* Calculates the the gradient for particle k */
    dwavefunction = zeros<rowvec>(nDimensions);
    rPlus = rMinus = r;
    wfCurrent = wavefunction(r);
    for (int i = 0; i < nDimensions; i++)
    {
        rPlus(k,i) += h;
        rMinus(k,i) -= h;
        wfMinus = wavefunction(rMinus);
        wfPlus = wavefunction(rPlus);
        dwavefunction(i) = wfPlus - wfMinus;
        rPlus(k,i) = rMinus(k,i) = r(k,i);
    }
    dwavefunction /= (2.0*wfCurrent*h);
//    dwavefunction /= (2.0*h); // not local

    return dwavefunction;
}

double Jastrow::localLaplacian(const int &k)
{
    /* Calculates the laplacian for particle k */
    lapl = 0.0;
    for (int i = 0; i < k; i++)
    {
        rki = rijNew(i,k);
        brki = 1.0 + beta*rki;
        lapl += a(k,i)/(rki*pow(brki, 3));
    }
    for (int i = k+1; i < nParticles; i++)
    {
        rki = rijNew(k,i);
        brki = 1.0 + beta*rki;
        lapl += a(k,i)/(rki*pow(brki, 3));
    }
    lapl *= 2.0;
    vec temp = localGradient(k);
    lapl += dot(temp, temp);

    return lapl;
}

double Jastrow::betaGradient(const int &k)
{
    /* Calculates the gradient for particle k */
    double delta = 0.0;
    for (int i = 0; i < k; i++)
    {
        rij = rijNew(i,k);
        delta -= a(i,k)*rij*rij/((beta*rij + 1.0)*(beta*rij + 1.0));
    }
    for (int i = k+1; i < nParticles; i++)
    {
        rij = rijNew(k,i);
        delta -= a(k,i)*rij*rij/((beta*rij + 1.0)*(beta*rij + 1.0));
    }

    return delta;
}

double Jastrow::localLaplacianNumerical(const mat &r, const int &k, const double &h)
{
    /* Calculates the laplacian for particle k */
    ddwavefunction = 0.0;
    rPlus = rMinus = r;
    wfCurrent = wavefunction(r);
    for (int i = 0; i < nDimensions; i++) {
        rPlus(k,i) += h;
        rMinus(k,i) -= h;
        wfMinus = wavefunction(rMinus);
        wfPlus = wavefunction(rPlus);
        ddwavefunction += wfMinus + wfPlus - 2.0*wfCurrent;
        rPlus(k,i) = rMinus(k,i) = r(k,i);
    }
    ddwavefunction /= (wfCurrent*h*h);

    return ddwavefunction;
}

void Jastrow::calculate_rij()
{
    rijOld.zeros();
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = i+1; j < nParticles; j++)
        {
            for (int k = 0; k < nDimensions; k++)
            {
                rijOld(i,j) += (rOld(i,k) - rOld(j,k))*(rOld(i,k) - rOld(j,k));
            }
            rijOld(i,j) = sqrt(rijOld(i,j));
        }
    }
}

void Jastrow::update_rij()
{
    int k = currentParticle;
    for (int i = 0; i < k; i++)
    {
        rijNew(i,k) = 0.0;
        for (int j = 0; j < nDimensions; j++)
        {
            rijNew(i,k) += (rNew(i,j) - rNew(k,j))*(rNew(i,j) - rNew(k,j));
        }
        rijNew(i,k) = sqrt(rijNew(i,k));
    }
    for (int i = k+1; i < nParticles; i++)
    {
        rijNew(k,i) = 0.0;
        for (int j = 0; j < nDimensions; j++)
        {
            rijNew(k,i) += (rNew(i,j) - rNew(k,j))*(rNew(i,j) - rNew(k,j));
        }
        rijNew(k,i) = sqrt(rijNew(k,i));
    }
}

void Jastrow::calculate_fij()
{
    double rij;
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < nParticles; j++)
        {
            rij = rijOld(i,j);
            fijOld(i,j) = a(i,j)*rij/(1.0 + beta*rij);
        }
    }
}

void Jastrow::update_fij()
{
    int k = currentParticle;
    double rij;
    for (int i = 0; i < k; i++)
    {
        rij = rijNew(i,k);
        fijNew(i,k) = a(i,k)*rij/(1.0 + beta*rij);
    }
    for (int i = k+1; i < nParticles; i++)
    {
        rij = rijNew(k,i);
        fijNew(k,i) = a(k,i)*rij/(1.0 + beta*rij);
    }
}


