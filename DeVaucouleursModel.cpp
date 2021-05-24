/*//////////////////////////////////////////////////////////////////
////     SpheCow -- Flexible dynamical models for galaxies      ////
////                and dark matter haloes                      ////
////     © Sterrenkundig Observatorium, Universiteit Gent       ////
///////////////////////////////////////////////////////////////// */

#include "DeVaucouleursModel.hpp"

//////////////////////////////////////////////////////////////////////

DeVaucouleursModel::DeVaucouleursModel(double Mtot, double Reff, const GaussLegendre* gl)
{
    _Mtot = Mtot;
    _Reff = Reff;
    _b = 7.669249442500804;
    _Sigmaff = pow(_b,8.0) / (40320.0*M_PI);
    _gl = gl;
}

//////////////////////////////////////////////////////////////////////

double DeVaucouleursModel::scale_radius() const
{
    return _Reff;
}

//////////////////////////////////////////////////////////////////////

double DeVaucouleursModel::surface_density(double R) const
{
    double dimf = _Mtot/pow(_Reff,2);
    double t = R/_Reff;
    double z = pow(t,0.25);
    return dimf * _Sigmaff * exp(-_b*z);
}

//////////////////////////////////////////////////////////////////////

double DeVaucouleursModel::derivative_surface_density(double R) const
{
    double dimf = _Mtot/pow(_Reff,3);
    double t = R/_Reff;
    double z = pow(t,0.25);
    double v1 = exp(-_b*z);
    double v2 = _b / 4.0 * (z/t);
    return -dimf * _Sigmaff * v1 * v2;
}

//////////////////////////////////////////////////////////////////////

double DeVaucouleursModel::second_derivative_surface_density(double R) const
{
    double dimf = _Mtot/pow(_Reff,4);
    double t = R/_Reff;
    double t2 = t*t;
    double z = pow(t,0.25);
    double v1 = exp(-_b*z);
    double v2 = _b * (3.0+_b*z) / 16.0 * (z/t2);
    return dimf * _Sigmaff * v1 * v2;
}

//////////////////////////////////////////////////////////////////////

double DeVaucouleursModel::third_derivative_surface_density(double R) const
{
    double dimf = _Mtot/pow(_Reff,5);
    double t = R/_Reff;
    double t3 = t*t*t;
    double z = pow(t,0.25);
    double v1 = exp(-_b*z);
    double v2 = _b * (21.0+9.0*_b*z+_b*_b*z*z) / 64.0 * (z/t3);
    return -dimf * _Sigmaff * v1 * v2;
}

//////////////////////////////////////////////////////////////////////
