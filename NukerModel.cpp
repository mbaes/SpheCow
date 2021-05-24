/*//////////////////////////////////////////////////////////////////
////     SpheCow -- Flexible dynamical models for galaxies      ////
////                and dark matter haloes                      ////
////     © Sterrenkundig Observatorium, Universiteit Gent       ////
///////////////////////////////////////////////////////////////// */

#include "NukerModel.hpp"

//////////////////////////////////////////////////////////////////////

NukerModel::NukerModel(double Mtot, double Rb, double alpha, double beta, double gamma, const GaussLegendre* gl)
{
    _Mtot = Mtot;
    _Rb = Rb;
    _alpha = alpha;
    _beta = beta;
    _gamma = gamma;
    double lg1 = lgamma((_beta-_gamma)/_alpha);
    double lg2 = lgamma((_beta-2.0)/_alpha);
    double lg3 = lgamma((2.0-_gamma)/_alpha);
    _Sigmab = _Mtot/(_Rb*_Rb) * exp(lg1-lg2-lg3) * _alpha / pow(2.0,(_beta-_gamma)/_alpha) / (2.0*M_PI);
    _gl = gl;
}

//////////////////////////////////////////////////////////////////////

double NukerModel::scale_radius() const
{
    return _Rb;
}

//////////////////////////////////////////////////////////////////////

double NukerModel::surface_density(double R) const
{
    double t = R/_Rb;
    double z = pow(t,_alpha);
    double q = (_beta-_gamma)/_alpha;
    return _Sigmab * pow(t,-_gamma) * pow((1.0+z)/2.0,-q);
}

//////////////////////////////////////////////////////////////////////

double NukerModel::derivative_surface_density(double R) const
{
    double t = R/_Rb;
    double z = pow(t,_alpha);
    double q = (_beta-_gamma)/_alpha;
    double ff = -pow(2.0,q) * _Sigmab/_Rb;
    return ff * pow(t,-1.0-_gamma) * pow(1.0+z,-1.0-q) * (_beta*z+_gamma);
}

//////////////////////////////////////////////////////////////////////

double NukerModel::second_derivative_surface_density(double R) const
{
    double t = R/_Rb;
    double z = pow(t,_alpha);
    double q = (_beta-_gamma)/_alpha;
    double ff = pow(2.0,q) * _Sigmab/(_Rb*_Rb);
    double v1 = pow(t,-2.0-_gamma);
    double v2 = pow(1.0+z,-2.0-q);
    double v3 = z*z*_beta*(1.0+_beta)
    + z*(_beta-_alpha*_beta+_gamma+_alpha*_gamma+2.0*_beta*_gamma)
    + _gamma*(1.0+_gamma);
    return ff * v1 * v2 * v3;
}

//////////////////////////////////////////////////////////////////////

double NukerModel::third_derivative_surface_density(double R) const
{
    double t = R/_Rb;
    double z = pow(t,_alpha);
    double ff = -pow(2.0,(_beta-_gamma)/_alpha) * _Sigmab/pow(_Rb,3);
    double v1 = pow(t,-3.0-_gamma);
    double v2 = pow(1.0+z,-3.0-(_beta-_gamma)/_alpha);
    double v3a = z*z*z*_beta*(1.0+_beta)*(2.0+_beta);
    double v3b = z*z*( _beta*(1.0-_alpha)*(4.0+_alpha+3.0*_beta)
                      +_gamma*(2.0+_alpha*_alpha+3.0*_alpha*(1.0+_beta)+3.0*_beta*(2.0+_beta)));
    double v3c = z*(-(1.0+_alpha)*(-4.0+_alpha-3.0*_gamma)*_gamma +
                    _beta*(2.0+_alpha*_alpha-3.0*_alpha*(1.0+_gamma)+3.0*_gamma*(2.0+_gamma)));
    double v3d = _gamma*(1.0+_gamma)*(2.0+_gamma);
    double v3 = v3a + v3b + v3c + v3d;
    return ff * v1 * v2 * v3;
}

//////////////////////////////////////////////////////////////////////
