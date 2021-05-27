/*//////////////////////////////////////////////////////////////////
////     SpheCow -- Flexible dynamical models for galaxies      ////
////                and dark matter haloes                      ////
////     © Sterrenkundig Observatorium, Universiteit Gent       ////
///////////////////////////////////////////////////////////////// */

#include "BurkertModel.hpp"

//////////////////////////////////////////////////////////////////////

BurkertModel::BurkertModel(double rhos, double rs, const GaussLegendre* gl)
{
    _rhos = rhos;
    _rs = rs;
    _gl = gl;
}

//////////////////////////////////////////////////////////////////////

double BurkertModel::scale_radius() const
{
    return _rhos;
}

//////////////////////////////////////////////////////////////////////

double BurkertModel::density(double r) const
{
    double dimf = _rhos;
    double t = r/_rs;
    double t2 = t*t;
    double z = (1.0+t) * (1.0+t2);
    return dimf / z;
}

//////////////////////////////////////////////////////////////////////

double BurkertModel::derivative_density(double r) const
{
    double dimf = _rhos/_rs;
    double t = r/_rs;
    double t2 = t*t;
    double z = (1.0+t) * (1.0+t2);
    return -dimf * (1.0 + 2.0*t + 3.0*t2) / (z*z);
}

//////////////////////////////////////////////////////////////////////

double BurkertModel::second_derivative_density(double r) const
{
    double dimf = _rhos/(_rs*_rs);
    double t = r/_rs;
    double t2 = t*t;
    double z = (1.0+t) * (1.0+t2);
    return dimf * 4.0*t2 * (3.0 + 4.0*t + 3.0*t2) / pow(z,3);
}

//////////////////////////////////////////////////////////////////////

double BurkertModel::mass(double r) const
{
    double dimf = _rhos*pow(_rs,3);
    double t = r/_rs;
    double t2 = t*t;
    return dimf * M_PI * (2.0*log(1.0+t) + log(1.0+t2) - 2.0*atan(t));
}

//////////////////////////////////////////////////////////////////////

double BurkertModel::potential(double r) const
{
    double dimf = _rhos*(_rs*_rs);
    double t = r/_rs;
    double t2 = t*t;
    double u = (1.0+t)/t;
    double v = (1.0-t)/t;
    return dimf * M_PI * (M_PI - 2.0*u*atan(t) + 2.0*u*log(1.0+t) + v*log(1.0+t2));
}

//////////////////////////////////////////////////////////////////////

