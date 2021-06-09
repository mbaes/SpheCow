/*//////////////////////////////////////////////////////////////////
////     SpheCow -- Flexible dynamical models for galaxies      ////
////                and dark matter haloes                      ////
////     © Sterrenkundig Observatorium, Universiteit Gent       ////
///////////////////////////////////////////////////////////////// */

#include "HypervirialModel.hpp"

//////////////////////////////////////////////////////////////////////

HypervirialModel::HypervirialModel(double Mtot, double rs, double p, const GaussLegendre* gl)
{
    _Mtot = Mtot;
    _rs = rs;
    _p = p;
    _gl = gl;
}

//////////////////////////////////////////////////////////////////////

double HypervirialModel::scale_radius() const
{
    return _rs;
}

//////////////////////////////////////////////////////////////////////

double HypervirialModel::density(double r) const
{
    double dimf = _Mtot/pow(_rs,3);
    double t = r/_rs;
    double tp = pow(t,_p);
    double v1 = pow(t,_p-2.0);
    double v2 = pow(1.0+tp,-2.0-1.0/_p);
    double v = v1 * v2;
    return dimf * (_p+1.0)/(4.0*M_PI) * v;
}

//////////////////////////////////////////////////////////////////////

double HypervirialModel::derivative_density(double r) const
{
    double dimf = _Mtot/pow(_rs,4);
    double t = r/_rs;
    double tp = pow(t,_p);
    double v1 = pow(t,_p-3.0);
    double v2 = pow(1.0+tp,-3.0-1.0/_p);
    double v3 = (2.0-_p+(3.0+_p)*tp);
    double v = v1 * v2 * v3;
    return -dimf * (_p+1.0)/(4.0*M_PI) * v;
}

//////////////////////////////////////////////////////////////////////

double HypervirialModel::second_derivative_density(double r) const
{
    double dimf = _Mtot/pow(_rs,5);
    double t = r/_rs;
    double tp = pow(t,_p);
    double v1 = pow(t,_p-4.0);
    double v2 = pow(1.0+tp,-4.0-1.0/_p);
    double p2 = _p*_p;
    double v3 = (6.0-5.0*_p+p2) + (17.0-3.0*_p-4.0*p2)*tp + (12.0+7.0*_p+p2)*tp*tp;
    double v = v1 * v2 * v3;
    return dimf * (_p+1.0)/(4.0*M_PI) * v;
}

//////////////////////////////////////////////////////////////////////

double HypervirialModel::mass(double r) const
{
    double dimf = _Mtot;
    double t = r/_rs;
    double q = 1.0+pow(t,-_p);
    return dimf / pow(q,1.0+1.0/_p);
}

//////////////////////////////////////////////////////////////////////

double HypervirialModel::total_mass() const
{
    return _Mtot;
}

//////////////////////////////////////////////////////////////////////

double HypervirialModel::potential(double r) const
{
    double dimf = _Mtot/_rs;
    double t = r/_rs;
    double z = 1.0+pow(t,_p);
    return dimf / pow(z,1.0/_p);
}

//////////////////////////////////////////////////////////////////////

double HypervirialModel::central_potential() const
{
    return _Mtot/_rs;
}

//////////////////////////////////////////////////////////////////////
