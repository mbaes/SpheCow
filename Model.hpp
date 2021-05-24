/*//////////////////////////////////////////////////////////////////
////     SpheCow -- Flexible dynamical models for galaxies      ////
////                and dark matter haloes                      ////
////     © Sterrenkundig Observatorium, Universiteit Gent       ////
///////////////////////////////////////////////////////////////// */

#ifndef MODEL_HPP
#define MODEL_HPP

#include "Basics.hpp"
class GaussLegendre;

//////////////////////////////////////////////////////////////////////

/** Model is the abstract base class for all spherical models. */

class Model
{
public:
 
    /** Virtual destructor of the Model class. */
    virtual ~Model() {};
    
    /** This pure virtual function returns the scale radius of the model. It is used in the numerical integration routines. */
    virtual double scale_radius() const = 0;
    
    /** This pure virtual function returns the total mass \f$M_{\text{tot}}\f$ . */
    virtual double total_mass() const = 0;
    
    /** This function returns the total potential energy \f$W_{\text{tot}}\f$ . It is calculated as \f[ W_{\text{tot}} = -4\pi G \int_0^\infty \rho(u)\, M(u)\, u\, {\text{d}} u. \f] The integration is performed using Gauss-Legendre quadrature. */
    double total_potential_energy() const;

    /** This pure virtual function returns the density \f$\rho(r)\f$ at radius \f$r\f$. */
    virtual double density(double r) const = 0;

    /** This pure virtual function returns the derivative of the density \f$\rho'(r)\f$ at radius \f$r\f$. */
    virtual double derivative_density(double r) const = 0;

    /** This pure virtual function returns the second derivative of the density \f$\rho''(r)\f$ at radius \f$r\f$. */
    virtual double second_derivative_density(double r) const = 0;

    /** This function returns the density slope \f$\gamma(r)\f$ at radius \f$r\f$. It is calculated as \f[ \gamma(r) = -\frac{{\text{d}}\log\rho}{{\text{d}}\log r}(r) = -\frac{r\,\rho'(r)}{\rho(r)}.\f] */
    double density_slope(double r) const;
    
    /** This pure virtual function returns the mass \f$M(r)\f$ at radius \f$r\f$. */
    virtual double mass(double r) const = 0;

    /** This function returns the circular velocity \f$v_{\text{c}}(r)\f$ at radius \f$r\f$. It is calculated as \f[ v_{\text{c}}(r) = \sqrt{\frac{GM(r)}{r}}.\f]  */
    double circular_velocity(double r) const;
    
    /** This pure virtual function returns the potential \f$\Psi(r)\f$ at radius \f$r\f$. */
    virtual double potential(double r) const = 0;
    
    /** This function returns the potential difference \f$\Psi(r_1)-\Psi(r_2)\f$ corresponding to two radii \f$r_1\f$ and \f$r_2\f$, with \f$r_2>r_1\f$. If these two radii are sufficiently apart, i.e., if \f$\epsilon \equiv r_2-r_1 > 10^{-4}\,r_{\text{s}}\f$, with \f$r_{\text{s}}\f$ the model scale radius, the routine directly uses the difference of the potential evaluated at the two radii. If \f$\epsilon \leq 10^{-4}\,r_{\text{s}}\f$, it uses the first terms in the Taylor expansion, \f[ \Psi(r_1)-\Psi(r_2) = -\left[\Psi(r_1+\epsilon)-\Psi(r_1)\right] \approx \frac{GM(r_1)}{r_1^2}\,\epsilon + \left[2\pi G\,\rho(r_1)-\frac{GM(r_1)}{r_1^3}\right] \epsilon^2.\f] */
    double potential_difference(double r1, double r2) const;

    /** This pure virtual function returns the surface density \f$\Sigma(R)\f$ at projected radius \f$R\f$. */
    virtual double surface_density(double R) const = 0;

    /** This pure virtual function returns the derivative of the surface density \f$\Sigma'(R)\f$ at projected radius \f$R\f$. */
    virtual double derivative_surface_density(double R) const = 0;

    /** This function returns the surface density slope \f$\gamma_{\text{p}}(R)\f$ at projected radius \f$R\f$. It is calculated as \f[ \gamma_{\text{p}}(R) = -\frac{{\text{d}}\log\Sigma}{{\text{d}}\log R}(R) = -\frac{R\,\Sigma'(R)}{\Sigma(R)}.\f] */
    double surface_density_slope(double R) const;
    
    /** This function returns the surface mass \f$M_{\text{p}}(R)\f$ at projected radius \f$R\f$. It is calculated as \f[ M_{\text{p}}(R) = 2\pi \int_0^R \Sigma(u)\,u\, {\text{d}}u. \f] The integration is performed using Gauss-Legendre quadrature. */
    double surface_mass(double R) const;

    /** This function returns the velocity dispersion \f$\sigma^2_{\text{iso}}(r)\f$ at radius \f$r\f$ under the assumption of an isotropic orbital structure. It is calculated as \f[ \sigma^2_{\text{iso}}(r) = \frac{G}{\rho(r)} \int_r^\infty \frac{\rho(u)\,M(u)\,{\text{d}} u}{u^2}.\f] The integration is performed using Gauss-Legendre quadrature. */
    double isotropic_dispersion(double r) const;
    
    /** This function returns the projected velocity dispersion \f$\sigma^2_{\text{p,iso}}(R)\f$ at projected radius \f$R\f$ under the assumption of an isotropic orbital structure. It is calculated as \f[ \sigma_{{\text{p}},{\text{iso}}}^2(R) = \frac{2G}{\Sigma(R)} \int_R^\infty \frac{\rho(u)\,M(u) \sqrt{u^2-R^2}\,{\text{d}} u}{u^2}.\f] The integration is performed using Gauss-Legendre quadrature. */
    double isotropic_projected_dispersion(double R) const;
    
    /** This function returns the distribution function \f$f_{\text{iso}}({\cal{E}})\f$ at binding energy \f${\cal{E}}=\Psi(r)\f$ under the assumption of an isotropic orbital structure. It is calculated as \f[ f_{\text{iso}}(\Psi(r)) = \frac{1}{2\sqrt2\,\pi^2} \int_r^\infty \frac{\Delta(u)\,{\text{d}}u}{\sqrt{\Psi(r)-\Psi(u)}},\f] with \f$\Delta(r)\f$ a function defined as \f[ \Delta(r) = \frac{r^2}{GM(r)}\left[\rho''(r) + \rho'(r) \left(\frac{2}{r} - \frac{4\pi\,\rho(r)\,r^2}{M(r)}\right) \right]. \f] The integration is performed using Gauss-Legendre quadrature. */
    virtual double isotropic_distribution_function(double r) const;
    
    /** This function returns the density \f$\rho(r)\f$ at radius \f$r\f$ calculated from the distribution function under the assumption of an isotropic orbital structure. It is calculated as \f[ \rho(r) = 4\sqrt2\,\pi\,G \int_r^\infty \frac{f_{\text{iso}}(\Psi(u))\,M(u) \sqrt{\Psi(r)-\Psi(u)}\,{\text{d}} u}{u^2}. \f] The integration is performed using Gauss-Legendre quadrature. This function can be used to check the implementation of new subclasses of the Model base class.*/
    double density_from_isotropic_distribution_function(double r) const;

    /** This function returns the velocity dispersion \f$\sigma^2_{\text{iso}}(r)\f$ at radius \f$r\f$ calculated from the distribution function under the assumption of an isotropic orbital structure. It is calculated as \f[ \sigma^2_{\text{iso}}(r) = \frac{8\sqrt2\,\pi}{3}\,\frac{G}{\rho(r)} \int_r^\infty \frac{f_{\text{iso}}(\Psi(u))\,M(u) [\Psi(r)-\Psi(u)]^{3/2}\,{\text{d}} u}{u^2}.\f] The integration is performed using Gauss-Legendre quadrature. This function can be used to check the implementation of new subclasses of the Model base class.*/
    double dispersion_from_isotropic_distribution_function(double r) const;

    /** This function returns the density-of-states function \f$g_{\text{iso}}({\cal{E}})\f$ at binding energy \f${\cal{E}}=\Psi(r)\f$ under the assumption of an isotropic orbital structure. It is calculated as \f[ g_{\text{iso}}(\Psi(r)) = 16\sqrt2\,\pi^2 \int_0^r u^2 \sqrt{\Psi(u)-\Psi(r)}\,{\text{d}} u.\f] The integration is performed using Gauss-Legendre quadrature. */
    double isotropic_density_of_states(double r) const;

    /** This function returns the total mass \f$M_{\text{tot}}\f$ calculated from the differential energy distribution \f${\cal{N}}({\cal{E}})\f$ under the assumption of an isotropic orbital structure. It is calculated as \f[ M_{\text{tot}} = G\int_0^\infty \frac{f_{\text{iso}}(\Psi(u))\, g_{\text{iso}}(\Psi(u))\, M(u)\,{\text{d}} u}{u^2}.\f] The integration is performed using Gauss-Legendre quadrature. This function can be used to check the implementation of new subclasses of the Model base class.*/
    double total_mass_from_isotropic_differential_energy_distribution() const;

    /** This function returns the total integrated binding energy \f${\cal{E}}_{\text{tot}}\f$ under the assumption of an isotropic orbital structure. It is calculated as \f[ {\cal{E}}_{\text{tot}} = \int_0^\infty \frac{ f_{\text{iso}}(\Psi(u))\, g_{\text{iso}}(\Psi(u))\, M(u)\, \Psi(u)\,{\text{d}} u}{u^2}. \f] The integration is performed using Gauss-Legendre quadrature.  */
    double isotropic_total_integrated_binding_energy() const;
    
    /** This function returns the total kinetic energy \f$T_{\text{tot}}\f$ under the assumption of an isotropic orbital structure. It is calculated as \f[ T_{\text{tot}} = 6\pi \int_0^\infty \rho(u)\,\sigma^2_{\text{iso}}(u)\,u^2\,{\text{d}} u,\f] with \f$\sigma^2_{\text{iso}}(r)\f$ the velocity dispersion. The integration is performed using Gauss-Legendre quadrature. */
    double isotropic_total_kinetic_energy() const;

    /** This function returns the radial velocity dispersion \f$\sigma^2_{r,\text{om}}(r)\f$ at radius \f$r\f$ under the assumption of an Osipkov-Merritt orbital structure with anisotropy radius \f$r_{\text{a}}\f$. It is calculated as \f[ \sigma_{r,\text{om}}^2(r) = \frac{G}{\rho(r)}\, \int_r^\infty \left(\frac{u^2+r_{\text{a}}^2}{r^2+r_{\text{a}}^2}\right) \frac{\rho(u)\,M(u)\,{\text{d}}u}{u^2}. \f] The integration is performed using Gauss-Legendre quadrature. */
    double osipkov_merritt_radial_dispersion(double r, double ra) const;
    
    /** This function returns the tangential velocity dispersion \f$\sigma^2_{\theta,\text{om}}(r) = \sigma_{\phi,{\text{om}}}^2(r)\f$ at radius \f$r\f$ under the assumption of an Osipkov-Merritt orbital structure with anisotropy radius \f$r_{\text{a}}\f$. It is calculated as \f[ \sigma_{\theta,{\text{om}}}^2(r) = \sigma_{\phi,{\text{om}}}^2(r) = \left(\frac{r_{\text{a}}^2}{r^2+r_{\text{a}}^2}\right) \sigma_{r,{\text{om}}}^2(r), \f] with \f$\sigma^2_{r,\text{om}}(r)\f$  the radial velocity dispersion. */
    double osipkov_merritt_tangential_dispersion(double r, double ra) const;
    
    /** This function returns the projected velocity dispersion \f$\sigma^2_{\text{p,om}}(R)\f$ at projected radius \f$R\f$ under the assumption of an Osipkov-Merritt orbital structure with anisotropy radius \f$r_{\text{a}}\f$. It is calculated as \f[ \sigma_{{\text{p}},{\text{om}}}^2(R) = \frac{G}{\Sigma(R)} \int_R^\infty \frac{w(u,R)\,\rho(u)\,M(u)\,{\text{d}}u}{u^2},\f] with \f[ w(u,R) = \left(\frac{u^2+r_{\text{a}}^2}{R^2+r_{\text{a}}^2}\right) \left( \frac{R^2+2r_{\text{a}}^2}{\sqrt{R^2+r_{\text{a}}^2}}\, \arctan \sqrt{\frac{u^2-R^2}{R^2+r_{\text{a}}^2}} - \frac{R^2\sqrt{u^2-R^2}}{u^2+r_{\text{a}}^2} \right). \f] The integration is performed using Gauss-Legendre quadrature. */
    double osipkov_merritt_projected_dispersion(double R, double ra) const;
    
    /** This function returns the distribution function \f$f_{\text{om}}(Q)\f$ at pseudo-binding energy \f$Q=\Psi(r)\f$ under the assumption of an Osipkov-Merritt orbital structure with anisotropy radius \f$r_{\text{a}}\f$. It is calculated as \f[ f_{\text{om}}(\Psi(r)) = \frac{1}{2\sqrt2\,\pi^2} \int_r^\infty \frac{\Delta_Q(u)\,{\text{d}}u}{\sqrt{\Psi(r)-\Psi(u)}},\f] with \f$\Delta_Q(r)\f$ a function defined as \f[ \Delta_Q(r) = \frac{r^2}{GM(r)}\left[\rho_Q''(r) + \rho_Q'(r) \left(\frac{2}{r} - \frac{4\pi\,\rho_Q(r)\,r^2}{M(r)}\right) \right], \f] with \f[ \rho_Q(r) = \left(1+ \frac{r^2}{r_{\text{a}}^2}\right) \rho(r).\f] The integration is performed using Gauss-Legendre quadrature. */
    virtual double osipkov_merritt_distribution_function(double r, double ra) const;

    /** This function returns the density \f$\rho(r)\f$ at radius \f$r\f$ calculated from the distribution function under the assumption of an Osipkov-Merritt orbital structure with anisotropy radius \f$r_{\text{a}}\f$. It is calculated as \f[ \rho(r) = 4\sqrt2\,\pi\,G \left(1+\frac{r^2}{r_{\text{a}}^2}\right)^{-1} \int_r^\infty \frac{f_{\text{om}}(\Psi(u))\,M(u) \sqrt{\Psi(r)-\Psi(u)}\,{\text{d}} u}{u^2}. \f] The integration is performed using Gauss-Legendre quadrature. This function can be used to check the implementation of new subclasses of the Model base class.*/
    double density_from_osipkov_merritt_distribution_function(double r, double ra) const;

    /** This function returns the radial velocity dispersion \f$\sigma^2_{r,\text{om}}(r)\f$ at radius \f$r\f$ calculated from the distribution function under the assumption of an Osipkov-Merritt orbital structure with anisotropy radius \f$r_{\text{a}}\f$. It is calculated as \f[ \sigma^2_{r,\text{om}}(r) = \frac{8\sqrt2\,\pi}{3}\,\frac{G}{\rho(r)} \left(1+\frac{r^2}{r_{\text{a}}^2}\right)^{-1} \int_r^\infty \frac{f_{\text{om}}(\Psi(u))\,M(u) [\Psi(r)-\Psi(u)]^{3/2}\,{\text{d}} u}{u^2}.\f] The integration is performed using Gauss-Legendre quadrature. This function can be used to check the implementation of new subclasses of the Model base class.*/
    double radial_dispersion_from_osipkov_merritt_distribution_function(double r, double ra) const;
    
    /** This function returns the pseudo-density-of-states function \f$g_{\text{om}}(Q)\f$ at pseudo-binding energy \f$Q=\Psi(r)\f$ under the assumption of an Osipkov-Merritt orbital structure with anisotropy radius \f$r_{\text{a}}\f$. It is calculated as \f[ g_{\text{om}}(\Psi(r)) = 16\sqrt2\,\pi^2 \int_0^r u^2 \left(1+\frac{u^2}{r_{\text{a}}^2}\right)^{-1} \sqrt{\Psi(u)-\Psi(r)}\,{\text{d}} u.\f] The integration is performed using Gauss-Legendre quadrature. */
    double osipkov_merritt_pseudo_density_of_states(double r, double ra) const;
    
    /** This function returns the total mass \f$M_{\text{tot}}\f$ calculated from the pseudo-differential energy distribution \f${\cal{N}}(Q)\f$ under the assumption of an Osipkov-Merritt orbital structure with anisotropy radius \f$r_{\text{a}}\f$. It is calculated as \f[ M_{\text{tot}} = G\int_0^\infty \frac{f_{\text{om}}(\Psi(u))\,g_{\text{om}}(\Psi(u))\,M(u)\,{\text{d}} u}{u^2}.\f] The integration is performed using Gauss-Legendre quadrature. This function can be used to check the implementation of new subclasses of the Model base class.*/
    double total_mass_from_osipkov_merritt_pseudo_differential_energy_distribution(double ra) const;
    
    /** This function returns the total kinetic energy \f$T_{\text{tot}}\f$ under the assumption of an Osipkov-Merritt orbital structure with anisotropy radius \f$r_{\text{a}}\f$. It is calculated as \f[ T_{\text{tot}} = 2\pi \int_0^\infty \left(\frac{u^2+3\,r_{\text{a}}^2}{u^2+r_{\text{a}}^2}\right) \rho(u)\,\sigma^2_{r,\text{om}}(u)\,u^2\,{\text{d}} u,\f] with \f$\sigma^2_{r,\text{om}}(r)\f$ the radial velocity dispersion. The integration is performed using Gauss-Legendre quadrature. */
    double osipkov_merritt_total_kinetic_energy(double ra) const;

protected:
    const GaussLegendre* _gl;
};

//////////////////////////////////////////////////////////////////////

#endif
