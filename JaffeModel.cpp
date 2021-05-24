/*//////////////////////////////////////////////////////////////////
////     SpheCow -- Flexible dynamical models for galaxies      ////
////                and dark matter haloes                      ////
////     © Sterrenkundig Observatorium, Universiteit Gent       ////
///////////////////////////////////////////////////////////////// */

#include "JaffeModel.hpp"

//////////////////////////////////////////////////////////////////////

JaffeModel::JaffeModel(double Mtot, double b, const GaussLegendre* gl)
{
    _Mtot = Mtot;
    _b = b;
    _gl = gl;
}

//////////////////////////////////////////////////////////////////////

double JaffeModel::scale_radius() const
{
    return _b;
}

//////////////////////////////////////////////////////////////////////

double JaffeModel::density(double r) const
{
    double dimf = _Mtot/pow(_b,3);
    double t = r/_b;
    double z = 1.0+t;
    return dimf / (4.0*M_PI) / (t*t*z*z);
}

//////////////////////////////////////////////////////////////////////

double JaffeModel::derivative_density(double r) const
{
    double dimf = _Mtot/pow(_b,4);
    double t = r/_b;
    double z = 1.0+t;
    return -dimf / (2.0*M_PI) * (1.0+2.0*t) / pow(t*z,3);
}

//////////////////////////////////////////////////////////////////////

double JaffeModel::second_derivative_density(double r) const
{
    double dimf = _Mtot/pow(_b,5);
    double t = r/_b;
    double z = 1.0+t;
    return dimf / (2.0*M_PI) * (3.0+10.0*t+10.0*t*t) / pow(t*z,4);
}

//////////////////////////////////////////////////////////////////////

double JaffeModel::mass(double r) const
{
    double dimf = _Mtot;
    double t = r/_b;
    return dimf * t/(1.0+t);
}

//////////////////////////////////////////////////////////////////////

double JaffeModel::potential(double r) const
{
    double dimf = _Mtot/_b;
    double t = r/_b;
    return dimf * log(1.0+1.0/t);
}

//////////////////////////////////////////////////////////////////////
