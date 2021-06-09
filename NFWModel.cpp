/*//////////////////////////////////////////////////////////////////
////     SpheCow -- Flexible dynamical models for galaxies      ////
////                and dark matter haloes                      ////
////     © Sterrenkundig Observatorium, Universiteit Gent       ////
///////////////////////////////////////////////////////////////// */

#include "NFWModel.hpp"

//////////////////////////////////////////////////////////////////////

NFWModel::NFWModel(double Mvir, double rs, double c, const GaussLegendre* gl)
{
    _Mvir = Mvir;
    _rs = rs;
    _c = c;
    _rvir = _rs*_c;
    _rhoff = 1.0 / (4.0*M_PI*(log(1.0+_c)-_c/(1.0+_c)));
    _gl = gl;
}

//////////////////////////////////////////////////////////////////////

double NFWModel::scale_radius() const
{
    return _rs;
}

//////////////////////////////////////////////////////////////////////

double NFWModel::density(double r) const
{
    double dimf = _Mvir/pow(_rs,3);
    double t = r/_rs;
    double z = 1.0+t;
    return dimf * _rhoff / (t*z*z);
}

//////////////////////////////////////////////////////////////////////

double NFWModel::derivative_density(double r) const
{
    double dimf = _Mvir/pow(_rs,4);
    double t = r/_rs;
    double z = 1.0+t;
    return -dimf * _rhoff * (1+3.0*t) / (t*t*pow(z,3));
}

//////////////////////////////////////////////////////////////////////

double NFWModel::second_derivative_density(double r) const
{
    double dimf = _Mvir/pow(_rs,5);
    double t = r/_rs;
    double z = 1.0+t;
    return dimf * _rhoff * 2.0*(1+4.0*t+6.0*t*t) / (t*t*t*pow(z,4));
}

 //////////////////////////////////////////////////////////////////////

double NFWModel::mass(double r) const
{
    double dimf = _Mvir;
    double t = r/_rs;
    return dimf * _rhoff * 4.0*M_PI * (log(1.0+t)-t/(1.0+t));
}

//////////////////////////////////////////////////////////////////////

double NFWModel::total_mass() const
{
    return std::numeric_limits<double>::infinity();
}

//////////////////////////////////////////////////////////////////////

double NFWModel::potential(double r) const
{
    double dimf = _Mvir/_rs;
    double t = r/_rs;
    return dimf * _rhoff * 4.0*M_PI * log(1.0+t)/t;
}

//////////////////////////////////////////////////////////////////////

double NFWModel::central_potential() const
{
    return _Mvir/_rs * _rhoff * 4.0*M_PI;
}

//////////////////////////////////////////////////////////////////////
