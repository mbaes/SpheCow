/*//////////////////////////////////////////////////////////////////
////     SpheCow -- Flexible dynamical models for galaxies      ////
////                and dark matter haloes                      ////
////     © Sterrenkundig Observatorium, Universiteit Gent       ////
///////////////////////////////////////////////////////////////// */

#include "SurfaceDensityModel.hpp"
#include "GaussLegendre.hpp"
#include <functional>

//////////////////////////////////////////////////////////////////////

double SurfaceDensityModel::total_mass() const
{
    std::function<double(double)> integrand = [&](double u) -> double
    {
        return surface_density(u) * u;
    };
    return 2.0*M_PI*_gl->integrate_0_infty(integrand,scale_radius());
}

//////////////////////////////////////////////////////////////////////

double SurfaceDensityModel::density(double r) const
{
    std::function<double(double)> integrand = [&](double u) -> double
    {
        return derivative_surface_density(u) / sqrt((u-r)*(u+r));
    };
    return -M_1_PI * _gl->integrate_r_infty(integrand,r,scale_radius());
}

//////////////////////////////////////////////////////////////////////

double SurfaceDensityModel::derivative_density(double r) const
{
    std::function<double(double)> integrand = [&](double u) -> double
    {
        return second_derivative_surface_density(u) * u/sqrt((u-r)*(u+r));
    };
    return -M_1_PI * _gl->integrate_r_infty(integrand,r,scale_radius())/r;
}

//////////////////////////////////////////////////////////////////////

double SurfaceDensityModel::second_derivative_density(double r) const
{
    std::function<double(double)> integrand = [&](double u) -> double
    {
        return third_derivative_surface_density(u) * u*u/sqrt((u-r)*(u+r));
    };
    return -M_1_PI * _gl->integrate_r_infty(integrand,r,scale_radius())/(r*r);
}

//////////////////////////////////////////////////////////////////////

double SurfaceDensityModel::mass(double r) const
{
    std::function<double(double)> integrand1 = [&](double u) -> double
    {
        return derivative_surface_density(u) * u*u;
    };
    double ans1 = -M_PI * _gl->integrate_0_r(integrand1,r,scale_radius());
    std::function<double(double)> integrand2 = [&](double u) -> double
    {
        double t = sqrt((u-r)*(u+r));
        return derivative_surface_density(u) * (u*u*atan(r/t)-r*t);
    };
    double ans2 = -2.0 * _gl->integrate_r_infty(integrand2,r,scale_radius());
    return ans1 + ans2;
}

//////////////////////////////////////////////////////////////////////

double SurfaceDensityModel::potential(double r) const
{
    std::function<double(double)> integrand1 = [&](double u) -> double
    {
        return derivative_surface_density(u) * u*u;
    };
    double ans1 = -M_PI/r * _gl->integrate_0_r(integrand1,r,scale_radius());
    std::function<double(double)> integrand2 = [&](double u) -> double
    {
        double t = sqrt((u-r)*(u+r));
        return derivative_surface_density(u) * (u*u*atan(r/t)+r*t);
    };
    double ans2 = -2.0/r * _gl->integrate_r_infty(integrand2,r,scale_radius());
    return ans1 + ans2;
}

//////////////////////////////////////////////////////////////////////

double SurfaceDensityModel::central_potential() const
{
    std::function<double(double)> integrand = [&](double u) -> double
    {
        return derivative_surface_density(u) * u;
    };
    return -4.0 * _gl->integrate_0_infty(integrand,scale_radius());
}

//////////////////////////////////////////////////////////////////////
