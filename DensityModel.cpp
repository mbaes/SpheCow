/*//////////////////////////////////////////////////////////////////
////     SpheCow -- Flexible dynamical models for galaxies      ////
////                and dark matter haloes                      ////
////     © Sterrenkundig Observatorium, Universiteit Gent       ////
///////////////////////////////////////////////////////////////// */

#include "DensityModel.hpp"
#include "GaussLegendre.hpp"
#include <functional>

//////////////////////////////////////////////////////////////////////

double DensityModel::total_mass() const
{
    std::function<double(double)> integrand = [&](double u) -> double
    {
        return density(u) * (u*u);
    };
    return 4.0*M_PI*_gl->integrate_0_infty(integrand,scale_radius());
}

//////////////////////////////////////////////////////////////////////

double DensityModel::mass(double r) const
{
    std::function<double(double)> integrand = [&](double u) -> double
    {
        return density(u) * (u*u);
    };
    return 4.0*M_PI * _gl->integrate_0_r(integrand,r,scale_radius());
}

//////////////////////////////////////////////////////////////////////

double DensityModel::potential(double r) const
{
    std::function<double(double)> integrand = [&](double u) -> double
    {
        return density(u) * u;
    };
    return mass(r)/r + 4.0*M_PI*_gl->integrate_r_infty(integrand,r,scale_radius());
}

//////////////////////////////////////////////////////////////////////

double DensityModel::surface_density(double R) const
{
    std::function<double(double)> integrand = [&](double u) -> double
    {
        return density(u) * u / sqrt((u-R)*(u+R));
    };
    return 2.0 * _gl->integrate_r_infty(integrand,R,scale_radius());
}

//////////////////////////////////////////////////////////////////////

double DensityModel::derivative_surface_density(double R) const
{
    std::function<double(double)> integrand = [&](double u) -> double
    {
        return (density(u)+derivative_density(u)*u) * u / sqrt((u-R)*(u+R));
    };
    return 2.0/R * _gl->integrate_r_infty(integrand,R,scale_radius());
}

//////////////////////////////////////////////////////////////////////
