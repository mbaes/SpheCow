/*//////////////////////////////////////////////////////////////////
////     SpheCow -- Flexible dynamical models for galaxies      ////
////                and dark matter haloes                      ////
////     © Sterrenkundig Observatorium, Universiteit Gent       ////
///////////////////////////////////////////////////////////////// */

#include "SigmoidSurfaceDensityModel.hpp"

//////////////////////////////////////////////////////////////////////

SigmoidSurfaceDensityModel::SigmoidSurfaceDensityModel(double Mtot, double Rb, double alpha, double beta, double gamma, const GaussLegendre* gl)
{
    _Mtot = Mtot;
    _Rb = Rb;
    _alpha = alpha;
    _beta = beta;
    _gamma = gamma;
    _Sigmac = 1.0;
    _gl = gl;
    _Sigmac = _Mtot / total_mass();
}

//////////////////////////////////////////////////////////////////////

double SigmoidSurfaceDensityModel::scale_radius() const
{
    return _Rb;
}

//////////////////////////////////////////////////////////////////////

double SigmoidSurfaceDensityModel::surface_density(double R) const
{
    double t = R/_Rb;
    double l = _alpha*log(t);
    double l2 = l*l;
    double s = sqrt(1.0+l2);
    double e = exp(-(_beta-_gamma)/(2.0*_alpha)*s);
    double ff = _Sigmac;
    return ff * e / pow(t,0.5*(_beta+_gamma));
}

//////////////////////////////////////////////////////////////////////

double SigmoidSurfaceDensityModel::derivative_surface_density(double R) const
{
    double t = R/_Rb;
    double l = _alpha*log(t);
    double l2 = l*l;
    double s = sqrt(1.0+l2);
    double e = exp(-(_beta-_gamma)/(2.0*_alpha)*s);
    double ff = -_Sigmac/_Rb;
    double v = (_beta-_gamma)*l + (_beta+_gamma)*s;
    return ff * e / (2.0*s) / pow(t,(_beta+_gamma+2.0)/2.0) * v;
}

//////////////////////////////////////////////////////////////////////

double SigmoidSurfaceDensityModel::second_derivative_surface_density(double R) const
{
    double t = R/_Rb;
    double l = _alpha*log(t);
    double l2 = l*l;
    double s = sqrt(1.0+l2);
    double e = exp(-(_beta-_gamma)/(2.0*_alpha)*s);
    double ff = _Sigmac/(_Rb*_Rb);
    double v1 = 2.0*_alpha*(_gamma-_beta);
    double v2 = 2.0*(_beta-_gamma)*(_beta+_gamma+1.0)*l*(1.0+l2);
    double v3 = (_beta+_gamma)*(_beta+_gamma+2.0)*s;
    double v4 = 2.0*(_beta*(_beta+1.0)+_gamma*(_gamma+1.0))*l2*s;
    double v = v1+v2+v3+v4;
    return ff * e / (4.0*s*s*s) / pow(t,0.5*(_beta+_gamma)+2.0) * v;
}

//////////////////////////////////////////////////////////////////////

double SigmoidSurfaceDensityModel::third_derivative_surface_density(double R) const
{
    double t = R/_Rb;
    double l = _alpha*log(t);
    double l2 = l*l;
    double s = sqrt(1.0+l2);
    double e = exp(-(_beta-_gamma)/(2.0*_alpha)*s);
    double ff = -_Sigmac/(_Rb*_Rb*_Rb);
    double v1 = (_beta-_gamma) * (16.0+7.0*_beta*_beta+2.0*_beta*(12.0+5.0*_gamma)+_gamma*(24.0+7.0*_gamma)) * pow(l,3);
    double v2 = 4.0 * (_beta-_gamma) * (_beta*_beta + (1.0+_gamma)*(2.0+_gamma) + _beta*(3.0+_gamma)) * pow(l,5);
    double v3 = 4.0 * (2.0+_beta+_gamma) * (_beta*(_beta+1.0) + _gamma*(_gamma+1.0) - _beta*_gamma) * pow(l,4) * s;
    double v4 = (2.0+_beta+_gamma) * (6.0*_alpha*(_gamma-_beta) + (_beta+_gamma)*(4.0+_beta+_gamma)*s);
    double v5 = (2.0+_beta+_gamma) * l*l * (6.0*_alpha*(-_beta+_gamma) + (_beta * (8.0+5.0*_beta) - 2.0*(-4.0+_beta)*_gamma + 5.0*_gamma*_gamma)*s);
    double v6 = (_beta-_gamma) * l * (8.0 + 3.0*_beta*_beta + 6.0*_beta*(2.0+_gamma) + 3.0*_gamma*(4.0+_gamma) - 6.0*_alpha*(2.0*_alpha +(_beta-_gamma) * s ));
    double v = v1+v2+v3+v4+v5+v6;
    return ff * e / (8.0*pow(s,5)) / pow(t,0.5*(_beta+_gamma)+3.0) * v;
}

//////////////////////////////////////////////////////////////////////
