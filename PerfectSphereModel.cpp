/*//////////////////////////////////////////////////////////////////
////     SpheCow -- Flexible dynamical models for galaxies      ////
////                and dark matter haloes                      ////
////     © Sterrenkundig Observatorium, Universiteit Gent       ////
///////////////////////////////////////////////////////////////// */

#include "PerfectSphereModel.hpp"

//////////////////////////////////////////////////////////////////////

PerfectSphereModel::PerfectSphereModel(double Mtot, double c, const GaussLegendre* gl)
{
    _Mtot = Mtot;
    _c = c;
    _gl = gl;
}

//////////////////////////////////////////////////////////////////////

double PerfectSphereModel::scale_radius() const
{
    return _c;
}

//////////////////////////////////////////////////////////////////////

double PerfectSphereModel::density(double r) const
{
    double dimf = _Mtot/pow(_c,3);
    double t = r/_c;
    double z = 1.0+t*t;
    return dimf * 1.0/(M_PI*M_PI) / (z*z);
}

//////////////////////////////////////////////////////////////////////

double PerfectSphereModel::derivative_density(double r) const
{
    double dimf = _Mtot/pow(_c,4);
    double t = r/_c;
    double z = 1.0+t*t;
    return -dimf * 4.0/(M_PI*M_PI) * t / pow(z,3);
}

//////////////////////////////////////////////////////////////////////

double PerfectSphereModel::second_derivative_density(double r) const
{
    double dimf = _Mtot/pow(_c,5);
    double t = r/_c;
    double z = 1.0+t*t;
    return dimf * 4.0/(M_PI*M_PI) * (5.0*t*t-1.0) / pow(z,4);
}

 //////////////////////////////////////////////////////////////////////

double PerfectSphereModel::mass(double r) const
{
    double dimf = _Mtot;
    double t = r/_c;
    return dimf * 2.0/M_PI * (atan(t)-t/(1.0+t*t));
}

//////////////////////////////////////////////////////////////////////

double PerfectSphereModel::potential(double r) const
{
    double dimf = _Mtot/_c;
    double t = r/_c;
    return dimf * 2.0/M_PI * atan(t) / t;
}

//////////////////////////////////////////////////////////////////////
