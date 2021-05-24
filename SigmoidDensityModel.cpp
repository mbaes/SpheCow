/*//////////////////////////////////////////////////////////////////
////     SpheCow -- Flexible dynamical models for galaxies      ////
////                and dark matter haloes                      ////
////     © Sterrenkundig Observatorium, Universiteit Gent       ////
///////////////////////////////////////////////////////////////// */

#include "SigmoidDensityModel.hpp"

//////////////////////////////////////////////////////////////////////

SigmoidDensityModel::SigmoidDensityModel(double Mtot, double rb, double alpha, double beta, double gamma, const GaussLegendre* gl)
{
    _Mtot = Mtot;
    _rb = rb;
    _alpha = alpha;
    _beta = beta;
    _gamma = gamma;
    _rhoc = 1.0;
    _gl = gl;
    _rhoc = _Mtot / total_mass();
}

//////////////////////////////////////////////////////////////////////

double SigmoidDensityModel::scale_radius() const
{
    return _rb;
}

//////////////////////////////////////////////////////////////////////

double SigmoidDensityModel::density(double r) const
{
    double t = r/_rb;
    double l = _alpha*log(t);
    double l2 = l*l;
    double s = sqrt(1.0+l2);
    double e = exp(-(_beta-_gamma)/(2.0*_alpha)*s);
    double ff = _rhoc;
    return ff * e / pow(t,0.5*(_beta+_gamma));
}

//////////////////////////////////////////////////////////////////////

double SigmoidDensityModel::derivative_density(double r) const
{
    double t = r/_rb;
    double l = _alpha*log(t);
    double l2 = l*l;
    double s = sqrt(1.0+l2);
    double e = exp(-(_beta-_gamma)/(2.0*_alpha)*s);
    double ff = -_rhoc/_rb;
    double v = (_beta-_gamma)*l + (_beta+_gamma)*s;
    return ff * e / (2.0*s) / pow(t,(_beta+_gamma+2.0)/2.0) * v;
}

//////////////////////////////////////////////////////////////////////

double SigmoidDensityModel::second_derivative_density(double r) const
{
    double t = r/_rb;
    double l = _alpha*log(t);
    double l2 = l*l;
    double s = sqrt(1.0+l2);
    double e = exp(-(_beta-_gamma)/(2.0*_alpha)*s);
    double ff = _rhoc/(_rb*_rb);
    double v1 = 2.0*_alpha*(_gamma-_beta);
    double v2 = 2.0*(_beta-_gamma)*(_beta+_gamma+1.0)*l*(1.0+l2);
    double v3 = (_beta+_gamma)*(_beta+_gamma+2.0)*s;
    double v4 = 2.0*(_beta*(_beta+1.0)+_gamma*(_gamma+1.0))*l2*s;
    double v = v1+v2+v3+v4;
    return ff * e / (4.0*s*s*s) / pow(t,0.5*(_beta+_gamma)+2.0) * v;
}

//////////////////////////////////////////////////////////////////////
