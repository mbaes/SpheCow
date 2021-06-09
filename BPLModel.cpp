/*//////////////////////////////////////////////////////////////////
////     SpheCow -- Flexible dynamical models for galaxies      ////
////                and dark matter haloes                      ////
////     © Sterrenkundig Observatorium, Universiteit Gent       ////
///////////////////////////////////////////////////////////////// */

#include "BPLModel.hpp"

//////////////////////////////////////////////////////////////////////

BPLModel::BPLModel(double Mtot, double rb, double beta, double gamma, const GaussLegendre* gl)
{
    _Mtot = Mtot;
    _rb = rb;
    _beta = beta;
    _gamma = gamma;
    _rhoff = (_beta-3.0)*(3.0-_gamma)/(_beta-_gamma)/(4.0*M_PI);
    _gl = gl;
}

//////////////////////////////////////////////////////////////////////

double BPLModel::scale_radius() const
{
    return _rb;
}

//////////////////////////////////////////////////////////////////////

double BPLModel::density(double r) const
{
    double dimf = _Mtot/pow(_rb,3);
    double t = r/_rb;
    double eta = (t<=1.0) ? _gamma : _beta;
    return dimf * _rhoff * pow(t,-eta);
}

//////////////////////////////////////////////////////////////////////

double BPLModel::derivative_density(double r) const
{
    double dimf = _Mtot/pow(_rb,4);
    double t = r/_rb;
    double eta = (t<=1.0) ? _gamma : _beta;
    return -dimf * _rhoff * eta * pow(t,-eta-1.0);
}

//////////////////////////////////////////////////////////////////////

double BPLModel::second_derivative_density(double r) const
{
    double dimf = _Mtot/pow(_rb,5);
    double t = r/_rb;
    double eta = (t<=1.0) ? _gamma : _beta;
    return dimf * _rhoff * eta * (eta+1.0) * pow(t,-eta-2.0);
}

//////////////////////////////////////////////////////////////////////

double BPLModel::mass(double r) const
{
    double dimf = _Mtot;
    double t = r/_rb;
    if (t<=1.0)
        return dimf * _rhoff * (4.0*M_PI) / (3.0-_gamma) * pow(t,3.0-_gamma);
    else
    {
        if (fabs(_beta-3.0)<1e-5)
            return dimf * _rhoff * (4.0*M_PI) * (1.0/(3.0-_gamma)+log(t));
        else
            return dimf * _rhoff * (4.0*M_PI) / (_beta-3.0) * ((_beta-_gamma)/(3.0-_gamma)-pow(t,3.0-_beta));
    }
}

//////////////////////////////////////////////////////////////////////

double BPLModel::total_mass() const
{
    return _Mtot;
}

//////////////////////////////////////////////////////////////////////

double BPLModel::potential(double r) const
{
    double dimf = _Mtot/_rb;
    double t = r/_rb;
    if (t<=1.0)
    {
        if (fabs(_gamma-2.0)<1e-5)
            return dimf * _rhoff * (4.0*M_PI) * ((_beta-1.0)/(_beta-2.0)-log(t));
        else
            return dimf * _rhoff * (4.0*M_PI) / (2.0-_gamma) * ((_beta-_gamma)/(_beta-2.0)-pow(t,2.0-_gamma)/(3.0-_gamma));
    }
    else
    {
        if (fabs(_beta-3.0)<1e-5)
            return dimf * _rhoff * (4.0*M_PI) / t * ((4.0-_gamma)/(3.0-_gamma)+log(t));
        else
            return dimf * _rhoff * (4.0*M_PI) / (_beta-3.0) / t * ((_beta-_gamma)/(3.0-_gamma)-pow(t,3.0-_beta)/(_beta-2.0));
    }
}

//////////////////////////////////////////////////////////////////////

double BPLModel::central_potential() const
{
    if (_gamma>=2.0)
        return std::numeric_limits<double>::infinity();
    else
        return (_Mtot/_rb) * _rhoff * (4.0*M_PI) * (_beta-_gamma) / ((2.0-_gamma)*(_beta-2.0));
}

//////////////////////////////////////////////////////////////////////

double BPLModel::isotropic_distribution_function(double r) const
{
    if (r<=_rb)
    {
        double ff = -1.0/(_rb*_rb) * (_beta-_gamma)*(3.0-_gamma) / (8.0*M_SQRT2*M_PI*M_PI*M_PI);
        double jump = ff/sqrt(potential_difference(r,_rb));
        return DensityModel::isotropic_distribution_function(r) + jump;
    }
    else
        return DensityModel::isotropic_distribution_function(r);
}

//////////////////////////////////////////////////////////////////////

double BPLModel::osipkov_merritt_distribution_function(double r, double ra) const
{
    if (r<=_rb)
    {
        double s = _rb/ra;
        double ff = -1.0/(_rb*_rb) * (_beta-_gamma)*(3.0-_gamma) / (8.0*M_SQRT2*M_PI*M_PI*M_PI) * (1.0+s*s);
        double jump = ff/sqrt(potential_difference(r,_rb));
        return DensityModel::osipkov_merritt_distribution_function(r,ra) + jump;
    }
    else
        return DensityModel::osipkov_merritt_distribution_function(r,ra);
}

//////////////////////////////////////////////////////////////////////
