/*//////////////////////////////////////////////////////////////////
////     SpheCow -- Flexible dynamical models for galaxies      ////
////                and dark matter haloes                      ////
////     © Sterrenkundig Observatorium, Universiteit Gent       ////
///////////////////////////////////////////////////////////////// */

#include "IsochroneModel.hpp"

//////////////////////////////////////////////////////////////////////

IsochroneModel::IsochroneModel(double Mtot, double b, const GaussLegendre* gl)
{
    _Mtot = Mtot;
    _b = b;
    _gl = gl;
}

//////////////////////////////////////////////////////////////////////

double IsochroneModel::scale_radius() const
{
    return _b;
}

//////////////////////////////////////////////////////////////////////

double IsochroneModel::density(double r) const
{
    double dimf = _Mtot/pow(_b,3);
    double t = r/_b;
    double t2 = t*t;
    double u = sqrt(1.0+t2);
    double p1 = 1.0 + 2.0*u;
    double p2 = 4.0*M_PI * pow(u,3) * pow(1.0+u,2);
    return dimf * p1/p2;
}

//////////////////////////////////////////////////////////////////////

double IsochroneModel::derivative_density(double r) const
{
    double dimf = _Mtot/pow(_b,4);
    double t = r/_b;
    double t2 = t*t;
    double u = sqrt(1.0+t2);
    double p1 = t * (11.0 + 8.0*t2 + 9.0*u);
    double p2 = 4.0*M_PI * pow(u,5) * pow(1.0+u,3);
    return -dimf * p1/p2;
}

//////////////////////////////////////////////////////////////////////

double IsochroneModel::second_derivative_density(double r) const
{
    double dimf = _Mtot/pow(_b,5);
    double t = r/_b;
    double t2 = t*t;
    double u = sqrt(1.0+t2);
    double p1 = 5.0 * (-4.0+13.0*t2+14.0*t2*t2 + 4.0*u*(-1.0+4.0*t2+2.0*t2*t2));
    double p2 = 4.0*M_PI * pow(u,7) * pow(1.0+u,4);
    return dimf * p1/p2;
}

 //////////////////////////////////////////////////////////////////////

double IsochroneModel::mass(double r) const
{
    double dimf = _Mtot;
    double t = r/_b;
    double t2 = t*t;
    double u = sqrt(1.0+t2);
    return dimf * (-2.0+(2.0+t2)/u)/t;
}

//////////////////////////////////////////////////////////////////////

double IsochroneModel::potential(double r) const
{
    double dimf = _Mtot/_b;
    double t = r/_b;
    double t2 = t*t;
    double u = sqrt(1.0+t2);
    return dimf / (1.0+u);
}

//////////////////////////////////////////////////////////////////////
