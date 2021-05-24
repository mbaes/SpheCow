/*//////////////////////////////////////////////////////////////////
////     SpheCow -- Flexible dynamical models for galaxies      ////
////                and dark matter haloes                      ////
////     © Sterrenkundig Observatorium, Universiteit Gent       ////
///////////////////////////////////////////////////////////////// */

#include "GammaModel.hpp"

//////////////////////////////////////////////////////////////////////

GammaModel::GammaModel(double Mtot, double b, double gamma, const GaussLegendre* gl)
{
    _Mtot = Mtot;
    _b = b;
    _gamma = gamma;
    _rhob = _Mtot/pow(_b,3) * (3.0-gamma)/(4.0*M_PI);
    _gl = gl;
}

//////////////////////////////////////////////////////////////////////

double GammaModel::scale_radius() const
{
    return _b;
}

//////////////////////////////////////////////////////////////////////

double GammaModel::density(double r) const
{
    double t = r/_b;
    return _rhob / (pow(t,_gamma)*pow(1.0+t,4.0-_gamma));
}

//////////////////////////////////////////////////////////////////////

double GammaModel::derivative_density(double r) const
{
    double t = r/_b;
    return -(_rhob/_b) * (4.0*t+_gamma) / (pow(t,_gamma+1.0)*pow(1.0+t,5.0-_gamma));
}

//////////////////////////////////////////////////////////////////////

double GammaModel::second_derivative_density(double r) const
{
    double t = r/_b;
    return _rhob/(_b*_b) * (20.0*t*t+10.0*t*_gamma+_gamma*(1.0+_gamma)) / (pow(t,_gamma+2.0)*pow(1.0+t,6.0-_gamma));
}

//////////////////////////////////////////////////////////////////////

double GammaModel::mass(double r) const
{
    double t = r/_b;
    return _Mtot * pow(t/(1.0+t),3.0-_gamma);
}

//////////////////////////////////////////////////////////////////////

double GammaModel::potential(double r) const
{
    double dimf = _Mtot/_b;
    double t = r/_b;
    double eps = 1e-3;
    if (fabs(_gamma-2.0)>eps)
        return dimf * 1.0/(2.0-_gamma) * (1.0-pow(t/(1.0+t),2.0-_gamma));
    double w = 2.0-_gamma;
    double q = log(t/(1.0+t));
    return dimf * (-q - 0.5*w*q*q - w*w/6.0*q*q*q);
}

//////////////////////////////////////////////////////////////////////
