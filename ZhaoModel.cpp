/*//////////////////////////////////////////////////////////////////
////     SpheCow -- Flexible dynamical models for galaxies      ////
////                and dark matter haloes                      ////
////     © Sterrenkundig Observatorium, Universiteit Gent       ////
///////////////////////////////////////////////////////////////// */

#include "ZhaoModel.hpp"

//////////////////////////////////////////////////////////////////////

ZhaoModel::ZhaoModel(double Mtot, double rb, double alpha, double beta, double gamma, const GaussLegendre* gl)
{
    _Mtot = Mtot;
    _rb = rb;
    _alpha = alpha;
    _beta = beta;
    _gamma = gamma;
    double lg1 = lgamma((_beta-_gamma)/_alpha);
    double lg2 = lgamma((_beta-3.0)/_alpha);
    double lg3 = lgamma((3.0-_gamma)/_alpha);
    _rhoff = _alpha * exp(lg1-lg2-lg3) / (4.0*M_PI);
    _gl = gl;
}

//////////////////////////////////////////////////////////////////////

double ZhaoModel::scale_radius() const
{
    return _rb;
}

//////////////////////////////////////////////////////////////////////

double ZhaoModel::density(double r) const
{
    double dimf = _Mtot/pow(_rb,3);
    double t = r/_rb;
    double z = pow(t,_alpha);
    double q = (_beta-_gamma)/_alpha;
    return dimf * _rhoff * pow(t,-_gamma) * pow(1.0+z,-q);
}

//////////////////////////////////////////////////////////////////////

double ZhaoModel::derivative_density(double r) const
{
    double dimf = _Mtot/pow(_rb,4);
    double t = r/_rb;
    double z = pow(t,_alpha);
    double q = (_beta-_gamma)/_alpha;
    double v1 = pow(t,-_gamma-1.0);
    double v2 = pow(1.0+z,-q-1.0);
    double v3 = _beta*z+_gamma;
    return -dimf * _rhoff * v1 * v2 * v3;
}

//////////////////////////////////////////////////////////////////////

double ZhaoModel::second_derivative_density(double r) const
{
    double dimf = _Mtot/pow(_rb,5);
    double t = r/_rb;
    double z = pow(t,_alpha);
    double q = (_beta-_gamma)/_alpha;
    double v1 = pow(t,-_gamma-2.0);
    double v2 = pow(1.0+z,-q-2.0);
    double v3 = _gamma*(_gamma+1.0) + z*((1.0+_alpha)*_gamma + _beta*(1.0-_alpha+2.0*_gamma)) + z*z*_beta*(_beta+1.0);
    return dimf * _rhoff * v1 * v2 * v3;
}

//////////////////////////////////////////////////////////////////////

double ZhaoModel::total_mass() const
{
    return _Mtot;
}

//////////////////////////////////////////////////////////////////////

double ZhaoModel::central_potential() const
{
    if (_gamma>=2.0)
        return std::numeric_limits<double>::infinity();
    else
    {
        double lg1 = lgamma((_beta-2.0)/_alpha);
        double lg2 = lgamma((2.0-_gamma)/_alpha);
        double lg3 = lgamma((_beta-3.0)/_alpha);
        double lg4 = lgamma((3.0-_gamma)/_alpha);
        return _Mtot/_rb * exp(lg1+lg2-lg3-lg4);
    }
}

//////////////////////////////////////////////////////////////////////
