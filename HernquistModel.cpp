/*//////////////////////////////////////////////////////////////////
////     SpheCow -- Flexible dynamical models for galaxies      ////
////                and dark matter haloes                      ////
////     © Sterrenkundig Observatorium, Universiteit Gent       ////
///////////////////////////////////////////////////////////////// */

#include "HernquistModel.hpp"

//////////////////////////////////////////////////////////////////////

HernquistModel::HernquistModel(double Mtot, double b, const GaussLegendre* gl)
{
    _Mtot = Mtot;
    _b = b;
    _gl = gl;
}

//////////////////////////////////////////////////////////////////////

double HernquistModel::scale_radius() const
{
    return _b;
}

//////////////////////////////////////////////////////////////////////

double HernquistModel::density(double r) const
{
    double dimf = _Mtot/pow(_b,3);
    double t = r/_b;
    double z = 1.0+t;
    return dimf * (0.5/M_PI) / (t*z*z*z);
}

//////////////////////////////////////////////////////////////////////

double HernquistModel::derivative_density(double r) const
{
    double dimf = _Mtot/pow(_b,4);
    double t = r/_b;
    double z = 1.0+t;
    return -dimf * (0.5/M_PI) * (1.0+4.0*t)/(t*t*pow(z,4));
}

//////////////////////////////////////////////////////////////////////

double HernquistModel::second_derivative_density(double r) const
{
    double dimf = _Mtot/pow(_b,5);
    double t = r/_b;
    double z = 1.0+t;
    return dimf / M_PI * (1.0+5.0*t+10.0*t*t)/(t*t*t*pow(z,5));
}

//////////////////////////////////////////////////////////////////////

double HernquistModel::mass(double r) const
{
    double dimf = _Mtot;
    double t = r/_b;
    double z = 1.0+t;
    return dimf * pow(t/z,2);
}

//////////////////////////////////////////////////////////////////////

double HernquistModel::potential(double r) const
{
    double dimf = _Mtot/_b;
    double t = r/_b;
    return dimf / (1.0+t);
}

//////////////////////////////////////////////////////////////////////
