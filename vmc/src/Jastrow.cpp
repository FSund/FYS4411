#include "Jastrow.h"

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

void Jastrow::initalize(const mat &r)
{
    rNew = rOld = r;

    calculate_rij();
    calculate_fij();
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

double Jastrow::wavefunction()
{
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

double Jastrow::getRatio()
{
    double dU = 0.0;
    int k = currentParticle;
    for (int i = 0; i < k; i++)
        dU += fijNew(i, k) - fijOld(i, k);
    for (int i = k+1; i < nParticles; i++)
        dU += fijNew(k, i) - fijOld(k, i);

    return exp(dU);
}

mat Jastrow::gradient()
{

}

double Jastrow::laplacian()
{

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

void Jastrow::calculate_rij()
{
    rijNew.zeros();
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = i+1; j < nParticles; j++)
        {
            for (int k = 0; k < nDimensions; k++)
            {
                rijNew(i,j) += (rNew(i,k) - rNew(j,k))*(rNew(i,k) - rNew(j,k));
            }
            rijNew(i,j) = sqrt(rijNew(i,j));
        }
    }
    rijOld = rijNew;
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
            rij = rijNew(i,j);
            fijNew(i,j) = a(i,j)*rij/(1.0 + beta*rij);
        }
    }

//    a = ((i+j)%2 == 0) ? 0.25 : 0.5;
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
//    a = ((i+k)%2 == 0) ? 0.25 : 0.5;
//    a = ((k+i)%2 == 0) ? 0.25 : 0.5;
}


