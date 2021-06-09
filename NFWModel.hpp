/*//////////////////////////////////////////////////////////////////
////     SpheCow -- Flexible dynamical models for galaxies      ////
////                and dark matter haloes                      ////
////     © Sterrenkundig Observatorium, Universiteit Gent       ////
///////////////////////////////////////////////////////////////// */

#ifndef NFWMODEL_HPP
#define NFWMODEL_HPP

#include "DensityModel.hpp"

//////////////////////////////////////////////////////////////////////

/** NFWModel is a subclass of the DensityModel class and represents spherical models with a Navarro, Frenk & White (NFW) density profile, \f[ \rho(r) = \frac{g(c)}{4\pi}\, \frac{M_{\text{vir}}}{r_{\text{s}}^3} \left(\frac{r}{r_{\text{s}}}\right)^{-1} \left(1+\frac{r}{r_{\text{s}}}\right)^{-2}, \f] with \f[ g(c) = \frac{1}{\log(1+c)-c/(1+c)},\f] and the scale radius \f$r_{\text{s}}\f$ given by \f[ r_{\text{s}} = \frac{r_{\text{vir}}}{c}.\f] The free parameters are the virial mass \f$M_{\text{vir}}\f$, the virial radius \f$r_{\text{vir}}\f$, and the concentration \f$c\f$. For more information, see <a href="https://ui.adsabs.harvard.edu/abs/2001MNRAS.321..155L/abstract">Łokas & Mamon (2001)</a>. */

class NFWModel : public DensityModel
{
public:
    
    /** Constructor of the NFWModel class. */
    NFWModel(double Mvir, double rvir, double c, const GaussLegendre* gl);
    
    /** This function returns the scale radius \f$r_{\text{s}}\f$ of the NFW model. */
    double scale_radius() const;

    /** This function returns the density \f$\rho(r)\f$ of the NFW model at radius \f$r\f$. */
    double density(double r) const;

    /** This function returns the derivative of the density \f$\rho'(r)\f$ of the NFW model at radius \f$r\f$. */
    double derivative_density(double r) const;

    /** This function returns the second derivative of the density \f$\rho''(r)\f$ of the NFW model at radius \f$r\f$. */
    double second_derivative_density(double r) const;

    /** This function returns the mass \f$M(r)\f$ of the NFW model at radius \f$r\f$. */
    double mass(double r) const;

    /** This function returns the total mass \f$M_{\text{tot}}\f$ of the NFW model. */
    double total_mass() const;

    /** This function returns the potential \f$\Psi(r)\f$ of the NFW model at radius \f$r\f$. */
    double potential(double r) const;

    /** This function returns the central potential \f$\Psi_0\f$ of the NFW model. */
    double central_potential() const;

private:
    
    /** The virial mass \f$M_{\text{vir}}\f$. */
    double _Mvir;

    /** The virial radius \f$r_{\text{vir}}\f$. */
    double _rvir;

    /** The concentration \f$c\f$. */
    double _c;

    /** The scale radius \f$r_{\text{s}}\f$. */
    double _rs;
    
    /** A pre-factor for the density.\f$. */
    double _rhoff;
};

//////////////////////////////////////////////////////////////////////

#endif
