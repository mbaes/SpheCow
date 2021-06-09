/*//////////////////////////////////////////////////////////////////
////     SpheCow -- Flexible dynamical models for galaxies      ////
////                and dark matter haloes                      ////
////     © Sterrenkundig Observatorium, Universiteit Gent       ////
///////////////////////////////////////////////////////////////// */

#include "EinastoModel.hpp"
#include <fstream>
#include <iostream>

//////////////////////////////////////////////////////////////////////

EinastoModel::EinastoModel(double Mtot, double rh, double n, const GaussLegendre* gl)
{
    _Mtot = Mtot;
    _rh = rh;
    _n = n;

    // Determine the value of d. Seach it among the exact values in the file Einastod.txt. This file contains exact values for all n between 0.01 and 15 with a spacing of 0.001. If it is not there exit.
    
    std::ifstream file("Einastod.txt");
    bool dset = false;
    double nf, df;
    for (int j=0; j<14990; j++)
    {
        file >> nf >> df;
        if (fabs(_n-nf)<1e-9)
        {
            _d = df;
            dset = true;
            break;
        }
    }
    if (!dset)
    {
        std::cerr << "Attempting to set up an EinastoModel object with unknown value for n"
        << std::endl;
        exit(1);
    }
    
    _rho0 = _Mtot/pow(_rh,3) * pow(_d,3.0*_n)/(4.0*M_PI*_n*tgamma(3.0*_n));
    _gl = gl;
}

//////////////////////////////////////////////////////////////////////

double EinastoModel::scale_radius() const
{
    return _rh;
}

//////////////////////////////////////////////////////////////////////

double EinastoModel::density(double r) const
{
    double t = r/_rh;
    double z = pow(t,1.0/_n);
    return _rho0 * exp(-_d*z);
}

//////////////////////////////////////////////////////////////////////

double EinastoModel::derivative_density(double r) const
{
    double t = r/_rh;
    double z = pow(t,1.0/_n);
    return -_rho0*_d/_n * exp(-_d*z) * z / r;
}

//////////////////////////////////////////////////////////////////////

double EinastoModel::second_derivative_density(double r) const
{
    double t = r/_rh;
    double z = pow(t,1.0/_n);
    return _rho0*_d/(_n*_n) * exp(-_d*z) * (_n-1.0+_d*z) * z / (r*r);
}

//////////////////////////////////////////////////////////////////////

double EinastoModel::total_mass() const
{
    return _Mtot;
}

//////////////////////////////////////////////////////////////////////

double EinastoModel::central_potential() const
{
    return _Mtot/_rh * pow(_d,_n) * tgamma(2.0*_n) / tgamma(3.0*_n);
}

//////////////////////////////////////////////////////////////////////
