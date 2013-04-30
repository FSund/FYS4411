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

rowvec Jastrow::localGradient(const int &k)
{
//    mat grad = zeros<mat>(nParticles, nDimensions);
//    double rij;
//    for (int k = 0; k < nParticles; k++)
//    {
//        for (int i = 0; i < k; i++)
//        {
//            rij = rijNew(i,k);
//            grad.row(k) += ((rNew.row(k) - rNew.row(i))/rij)*
//                               ( a(i,k)/(1 + beta*rij)*(1 + beta*rij) );
//        }
//        for (int i = k+1; i < nParticles; i++)
//        {
//            rij = rijNew(k,i);
//            grad.row(k) -= ((rNew.row(i) - rNew.row(k))/rij)*
//                               ( a(k,i)/(1 + beta*rij)*(1 + beta*rij) );
//        }
//    }

    // lecture notes p. 515-517
    // only for particle k!
    grad = zeros<rowvec>(nDimensions);
    double rij;
    for (int i = 0; i < k; i++)
    {
        rij = rijNew(i,k);
        grad += (rNew.row(k) - rNew.row(i)) * a(i,k)/(rij*pow((1 + beta*rij),2));
    }
    for (int i = k+1; i < nParticles; i++)
    {
        rij = rijNew(k,i);
        grad -= (rNew.row(i) - rNew.row(k)) * a(k,i)/(rij*pow((1 + beta*rij),2));
    }

    return grad;
}

//mat Jastrow::localGradientNumerical(const double &h)
//{
//    mat gradient = zeros<rowvec>(nParticles, nDimensions);
//    double wfCurrent, wfMinus, wfPlus;

//    // computing the first derivative
//    rPlus = rMinus = rNew;
//    wfCurrent = wavefunction();
//    for (int i = 0; i < nParticles; i++) {
//        for (int j = 0; j < nDimensions; j++) {
//            rPlus(i,j) += h;
//            rMinus(i,j) -= h;
//            wfMinus = wavefunction(rMinus);
//            wfPlus = wavefunction(rPlus);
//            gradient(i,j) = wfPlus - wfMinus;
//            rPlus(i,j) = rMinus(i,j) = rNew(i,j);
//        }
//    }
//    gradient /= (2.0*wfCurrent*h);

//    return gradient;
//}

double Jastrow::localLaplacian(const int &k)
{
//    lapl = zeros<rowvec>(nParticles);
//    double arg, rki, rkj;
//    for (int k = 0; k < nParticles; k++)
//    {
//        for (int j = 0; j < nParticles; j++)
//        {
//            if (j == k) continue;
//            rkj = 0.0;
//            for (int l = 0; l < nDimensions; l++)
//                rkj += rNew(k,l) - rNew(k,j);
//            rkj = sqrt(rkj);
//            lapl(k) += 2.0*a(k,j)*(1.0/rkj + 1.0/(1.0 + beta*rkj))/pow((1 + beta*rkj),2);
//            for (int i = 0; i < nParticles; i++)
//            {
//                if (i == k) continue;
//                arg = rki = 0.0;
//                for (int l = 0; l < nDimensions; l++)
//                {
//                    arg += (rNew(k,l) - rNew(i,l))*(rNew(k,l) - rNew(j,l)); // dot product
//                    rki += rNew(k,l) - rNew(k,i);
//                }
//                rki = sqrt(rki);
//                lapl(k) += arg
//                    *(a(k,i)/(rki*pow((1.0 + beta*rki),2)))
//                    *(a(k,j)/(rkj*pow((1.0 + beta*rkj),2)));
//            }
//        }
//    }

//    lapl = 0.0;
//    double arg, rki, rkj;
//    for (int j = 0; j < nParticles; j++)
//    {
//        if (j == k) continue;
//        rkj = 0.0;
//        for (int l = 0; l < nDimensions; l++)
//            rkj += rNew(k,l) - rNew(k,j);
//        rkj = sqrt(rkj);
//        lapl += 2.0*a(k,j)*(1.0/rkj + 1.0/(1.0 + beta*rkj))/pow((1 + beta*rkj),2);
//        for (int i = 0; i < nParticles; i++)
//        {
//            if (i == k) continue;
//            arg = rki = 0.0;
//            for (int l = 0; l < nDimensions; l++)
//            {
//                arg += (rNew(k,l) - rNew(i,l))*(rNew(k,l) - rNew(j,l)); // dot product
//                rki += rNew(k,l) - rNew(k,i);
//            }
//            rki = sqrt(rki);
//            lapl += arg
//                    *(a(k,i)/(rki*pow((1.0 + beta*rki),2)))
//                    *(a(k,j)/(rkj*pow((1.0 + beta*rkj),2)));
//        }
//    }

//    return lapl;

    // https://github.com/sigvebs/VMC/blob/master/QD/QD_Jastrow.cpp
    lapl = 0.0;
    double rki;
    double brki;
    for (int i = 0; i < k; i++)
    {
        rki = 0.0;
        for (int l = 0; l < nDimensions; l++)
            rki += (rNew(k,l) - rNew(i,l))*(rNew(k,l) - rNew(i,l));
        rki = sqrt(rki);
        brki = 1.0 + beta*rki;
        lapl += a(k,i)*brki/(rki*pow(brki, 3));
    }
    for (int i = k+1; i < nParticles; i++)
    {
        rki = 0.0;
        for (int l = 0; l < nDimensions; l++)
            rki += (rNew(k,l) - rNew(i,l))*(rNew(k,l) - rNew(i,l));
        rki = sqrt(rki);
        brki = 1.0 + beta*rki;
        lapl += a(k,i)*brki/(rki*pow(brki, 3));
    }
    lapl += dot(localGradient(k), localGradient(k)); // inefficient!!
    // lapl += dot(grad, grad); // should have calculated the gradient before this

//    cout << "lapl jastrow = " << lapl << endl;
    return lapl;
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


