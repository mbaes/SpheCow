/*//////////////////////////////////////////////////////////////////
////     SpheCow -- Flexible dynamical models for galaxies      ////
////                and dark matter haloes                      ////
////     © Sterrenkundig Observatorium, Universiteit Gent       ////
///////////////////////////////////////////////////////////////// */

#ifndef BURKERTMODEL_HPP
#define BURKERTMODEL_HPP

#include "DensityModel.hpp"

//////////////////////////////////////////////////////////////////////

/** BurkertModel is a subclass of the DensityModel class and represents spherical models with a Burkert density profile, \f[ \rho(r) = \rho_{\text{s}} \left(\frac{r}{r_{\text{s}}}\right)^{-1} \left(1+ \frac{r^2}{r_{\text{s}}^2}\right)^{-1}. \f] The free parameters are the density scale \f$\rho_{\text{s}}\f$ and the scale length \f$r_{\text{s}}\f$. For more information, see <a href="https://ui.adsabs.harvard.edu/abs/1995ApJ...447L..25B/abstract">Burkert (1995)</a>. */

class BurkertModel : public DensityModel
{
public:

    /** Constructor of the BurkertModel class. */
    BurkertModel(double rhos, double rs, const GaussLegendre* gl);

    /** This function returns the scale radius \f$r_{\text{s}}\f$ of the Burkert model. */
    double scale_radius() const;

    /** This function returns the density \f$\rho(r)\f$ of the Burkert model at radius \f$r\f$. */
    double density(double r) const;

    /** This function returns the derivative of the density \f$\rho'(r)\f$ of the Burkert model at radius \f$r\f$. */
    double derivative_density(double r) const;

    /** This function returns the second derivative of the density \f$\rho''(r)\f$ of the Burkert model at radius \f$r\f$. */
    double second_derivative_density(double r) const;

    /** This function returns the mass \f$M(r)\f$ of the Burkert model at radius \f$r\f$. */
    double mass(double r) const;

    /** This function returns the total mass \f$M_{\text{tot}}\f$ of the Burkert model. */
    double total_mass() const;

    /** This function returns the potential \f$\Psi(r)\f$ of the Burkert model at radius \f$r\f$. */
    double potential(double r) const;

    /** This function returns the central potential \f$\Psi_0\f$ of the Burkert model. */
    double central_potential() const;

private:
    
    /** The density scale \f$\rho_{\text{s}}\f$. */
    double _rhos;

    /** The scale radius \f$r_{\text{s}}\f$. */
    double _rs;
};

//////////////////////////////////////////////////////////////////////

#endif
