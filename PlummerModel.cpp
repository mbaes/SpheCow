/*//////////////////////////////////////////////////////////////////
////     SpheCow -- Flexible dynamical models for galaxies      ////
////                and dark matter haloes                      ////
////     © Sterrenkundig Observatorium, Universiteit Gent       ////
///////////////////////////////////////////////////////////////// */

#include "PlummerModel.hpp"

//////////////////////////////////////////////////////////////////////

PlummerModel::PlummerModel(double Mtot, double c, const GaussLegendre* gl)
{
    _Mtot = Mtot;
    _c = c;
    _gl = gl;
}

//////////////////////////////////////////////////////////////////////

double PlummerModel::scale_radius() const
{
    return _c;
}

//////////////////////////////////////////////////////////////////////

double PlummerModel::density(double r) const
{
    double dimf = _Mtot/pow(_c,3);
    double t = r/_c;
    double z = sqrt(1.0+t*t);
    return dimf * 3.0/(4.0*M_PI) / pow(z,5);
}

//////////////////////////////////////////////////////////////////////

double PlummerModel::derivative_density(double r) const
{
    double dimf = _Mtot/pow(_c,4);
    double t = r/_c;
    double z = sqrt(1.0+t*t);
    return -dimf * 15.0/(4.0*M_PI) * t/pow(z,7);
}

//////////////////////////////////////////////////////////////////////

double PlummerModel::second_derivative_density(double r) const
{
    double dimf = _Mtot/pow(_c,5);
    double t = r/_c;
    double z = sqrt(1.0+t*t);
    return dimf * 15.0/(4.0*M_PI) * (6.0*t*t-1.0) / pow(z,9);
}

//////////////////////////////////////////////////////////////////////

double PlummerModel::mass(double r) const
{
    double dimf = _Mtot;
    double t = r/_c;
    double z = sqrt(1.0+t*t);
    return dimf * pow(t/z,3);
}

//////////////////////////////////////////////////////////////////////

double PlummerModel::total_mass() const
{
    return _Mtot;
}

//////////////////////////////////////////////////////////////////////

double PlummerModel::potential(double r) const
{
    double dimf = _Mtot/_c;
    double t = r/_c;
    return dimf / sqrt(1.0+t*t);
}

//////////////////////////////////////////////////////////////////////

double PlummerModel::central_potential() const
{
    return _Mtot/_c;
}

//////////////////////////////////////////////////////////////////////

