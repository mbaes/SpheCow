/*//////////////////////////////////////////////////////////////////
////     SpheCow -- Flexible dynamical models for galaxies      ////
////                and dark matter haloes                      ////
////     © Sterrenkundig Observatorium, Universiteit Gent       ////
///////////////////////////////////////////////////////////////// */

#include "Model.hpp"
#include "GaussLegendre.hpp"
#include <functional>

//////////////////////////////////////////////////////////////////////

double Model::total_potential_energy() const
{
    std::function<double(double)> integrand = [&](double u) -> double
    {
        return density(u) * mass(u) * u;
    };
    return -4.0*M_PI * _gl->integrate_0_infty(integrand,scale_radius());
}

//////////////////////////////////////////////////////////////////////

double Model::density_slope(double r) const
{
    return -r * derivative_density(r) / density(r);
}

//////////////////////////////////////////////////////////////////////

double Model::circular_velocity(double r) const
{
    return mass(r) / r;
}

//////////////////////////////////////////////////////////////////////

double Model::potential_difference(double r1, double r2) const
{
    double eps = r2 - r1;
    if (eps>1e-4*scale_radius()) return potential(r1)-potential(r2);
    double M = mass(r1);
    double rho = density(r1);
    return M/(r1*r1)*eps + (2.0*M_PI*rho-M/pow(r1,3))*eps*eps;
}

//////////////////////////////////////////////////////////////////////

double Model::surface_density_slope(double R) const
{
    return -R * derivative_surface_density(R) / surface_density(R);
}

//////////////////////////////////////////////////////////////////////

double Model::surface_mass(double R) const
{
    std::function<double(double)> integrand = [&](double u) -> double
    {
        return surface_density(u) * u;
    };
    return 2.0*M_PI*_gl->integrate_0_r(integrand,R,scale_radius());
}

//////////////////////////////////////////////////////////////////////

double Model::isotropic_dispersion(double r) const
{
    std::function<double(double)> integrand = [&](double u) -> double
    {
        return density(u) * mass(u) / (u*u);
    };
    return _gl->integrate_r_infty(integrand,r,scale_radius()) / density(r);
}

//////////////////////////////////////////////////////////////////////

double Model::isotropic_projected_dispersion(double R) const
{
    std::function<double(double)> integrand = [&](double u) -> double
    {
        return density(u) * mass(u) / (u*u) * sqrt((u-R)*(u+R));
    };
    return 2.0*_gl->integrate_r_infty(integrand,R,scale_radius()) / surface_density(R);
}

//////////////////////////////////////////////////////////////////////

double Model::isotropic_distribution_function(double r) const
{
    std::function<double(double)> integrand = [&](double u) -> double
    {
        double M = mass(u);
        double rho = density(u);
        double drho = derivative_density(u);
        double d2rho = second_derivative_density(u);
        double Delta = u*u/M * (d2rho + drho*(2.0/u-4.0*M_PI*rho*u*u/M));
        return Delta/sqrt(fabs(potential_difference(r,u)));
    };
    return 1.0/(2.0*M_SQRT2*M_PI*M_PI) * _gl->integrate_r_infty(integrand,r,scale_radius());
}

//////////////////////////////////////////////////////////////////////

double Model::density_from_isotropic_distribution_function(double r) const
{
    std::function<double(double)> integrand = [&](double u) -> double
    {
        double z = sqrt(fabs(potential_difference(r,u)));
        return isotropic_distribution_function(u) * mass(u) * z / (u*u);
    };
    return 4.0*M_SQRT2*M_PI * _gl->integrate_r_infty(integrand,r,scale_radius());
}

//////////////////////////////////////////////////////////////////////

double Model::dispersion_from_isotropic_distribution_function(double r) const
{
    std::function<double(double)> integrand = [&](double u) -> double
    {
        double z = sqrt(fabs(potential_difference(r,u)));
        return isotropic_distribution_function(u) * mass(u) * (z*z*z) / (u*u);
    };
    return 8.0*M_SQRT2*M_PI/3.0 * _gl->integrate_r_infty(integrand,r,scale_radius()) / density(r);
}

//////////////////////////////////////////////////////////////////////

double Model::isotropic_density_of_states(double r) const
{
    std::function<double(double)> integrand = [&](double u) -> double
    {
        return (u*u) * sqrt(fabs(potential_difference(u,r)));
    };
    return 16.0*M_SQRT2*M_PI*M_PI * _gl->integrate_0_r(integrand,r,scale_radius());
}

//////////////////////////////////////////////////////////////////////

double Model::total_mass_from_isotropic_differential_energy_distribution() const
{
    std::function<double(double)> integrand = [&](double u) -> double
    {
        double df = isotropic_distribution_function(u);
        double g = isotropic_density_of_states(u);
        return df * g * mass(u) / (u*u);
    };
    return _gl->integrate_0_infty(integrand,scale_radius());
}

//////////////////////////////////////////////////////////////////////

double Model::isotropic_total_integrated_binding_energy() const
{
    std::function<double(double)> integrand = [&](double u) -> double
    {
        double df = isotropic_distribution_function(u);
        double g = isotropic_density_of_states(u);
        return df * g * mass(u) * potential(u) / (u*u);
    };
    return _gl->integrate_0_infty(integrand,scale_radius());
}

//////////////////////////////////////////////////////////////////////

double Model::isotropic_total_kinetic_energy() const
{
    std::function<double(double)> integrand = [&](double u) -> double
    {
        return density(u) * isotropic_dispersion(u) * (u*u);
    };
    return 6.0*M_PI * _gl->integrate_0_infty(integrand,scale_radius());
}

//////////////////////////////////////////////////////////////////////

double Model::osipkov_merritt_radial_dispersion(double r, double ra) const
{
    std::function<double(double)> integrand = [&](double u) -> double
    {
        
        double rhoQ = density(u) * (1.0+u*u/(ra*ra));
        return rhoQ * mass(u) / (u*u);
    };
    return _gl->integrate_r_infty(integrand,r,scale_radius()) / (1.0+r*r/(ra*ra)) / density(r);
}

//////////////////////////////////////////////////////////////////////

double Model::osipkov_merritt_tangential_dispersion(double r, double ra) const
{
    return 1.0/(1.0+r*r/(ra*ra)) * osipkov_merritt_radial_dispersion(r,ra);
}

//////////////////////////////////////////////////////////////////////

double Model::osipkov_merritt_projected_dispersion(double R, double ra) const
{
    std::function<double(double)> integrand = [&](double u) -> double
    {
        double f = (u*u+ra*ra) / (R*R+ra*ra);
        double t1 = (R*R+2.0*ra*ra) / sqrt(R*R+ra*ra) * atan(sqrt((u-R)*(u+R)/(R*R+ra*ra)));
        double t2 = -R*R * sqrt((u-R)*(u+R))/(u*u+ra*ra);
        double w = f * (t1+t2);
        return w * density(u) * mass(u) / (u*u);
    };
    return _gl->integrate_r_infty(integrand,R,scale_radius()) / surface_density(R);
}

//////////////////////////////////////////////////////////////////////

double Model::osipkov_merritt_distribution_function(double r, double ra) const
{
    std::function<double(double)> integrand = [&](double u) -> double
    {
        double M = mass(u);
        double rho = density(u);
        double drho = derivative_density(u);
        double d2rho = second_derivative_density(u);
        double z = (1.0+u*u/(ra*ra));
        double drhoQ = 2.0 *u/(ra*ra)*rho + z*drho;
        double d2rhoQ = 2.0*rho/(ra*ra) + 4.0*u/(ra*ra)*drho + z*d2rho;
        double DeltaQ = u*u/M * (d2rhoQ + drhoQ*(2.0/u-4.0*M_PI*rho*u*u/M));
        return DeltaQ / sqrt(fabs(potential_difference(r,u)));
    };
    return 1.0/(2.0*M_SQRT2*M_PI*M_PI) * _gl->integrate_r_infty(integrand,r,scale_radius());
}

//////////////////////////////////////////////////////////////////////

double Model::density_from_osipkov_merritt_distribution_function(double r, double ra) const
{
    std::function<double(double)> integrand = [&](double u) -> double
    {
        double z = sqrt(fabs(potential_difference(r,u)));
        return osipkov_merritt_distribution_function(u,ra) * mass(u) * z / (u*u);
    };
    return 4.0*M_SQRT2*M_PI / (1.0+r*r/(ra*ra)) * _gl->integrate_r_infty(integrand,r,scale_radius());
}

//////////////////////////////////////////////////////////////////////

double Model::radial_dispersion_from_osipkov_merritt_distribution_function(double r, double ra) const
{
    std::function<double(double)> integrand = [&](double u) -> double
    {
        double z = sqrt(fabs(potential_difference(r,u)));
        return osipkov_merritt_distribution_function(u,ra) * mass(u) * (z*z*z) / (u*u);
    };
    return 8.0*M_SQRT2*M_PI/3.0 / (1.0+r*r/(ra*ra)) * _gl->integrate_r_infty(integrand,r,scale_radius()) / density(r);
}

//////////////////////////////////////////////////////////////////////

double Model::osipkov_merritt_pseudo_density_of_states(double r, double ra) const
{
    std::function<double(double)> integrand = [&](double u) -> double
    {
        return u*u/(1.0+u*u/(ra*ra)) * sqrt(fabs(potential_difference(u,r)));
    };
    return 16.0*M_SQRT2*M_PI*M_PI * _gl->integrate_0_r(integrand,r,scale_radius());
}

//////////////////////////////////////////////////////////////////////

double Model::total_mass_from_osipkov_merritt_pseudo_differential_energy_distribution(double ra) const
{
    std::function<double(double)> integrand = [&](double u) -> double
    {
        double df = osipkov_merritt_distribution_function(u,ra);
        double g = osipkov_merritt_pseudo_density_of_states(u,ra);
        return df * g * mass(u) / (u*u);
    };
    return _gl->integrate_0_infty(integrand,scale_radius());
}

//////////////////////////////////////////////////////////////////////

double Model::osipkov_merritt_total_kinetic_energy(double ra) const
{
    std::function<double(double)> integrand = [&](double u) -> double
    {
        return (1.0+2.0*ra*ra/(u*u+ra*ra)) * density(u) * osipkov_merritt_radial_dispersion(u,ra) * (u*u);
    };
    return 2.0*M_PI*_gl->integrate_0_infty(integrand,scale_radius());
}

//////////////////////////////////////////////////////////////////////
