/*//////////////////////////////////////////////////////////////////
////     SpheCow -- Flexible dynamical models for galaxies      ////
////                and dark matter haloes                      ////
////     © Sterrenkundig Observatorium, Universiteit Gent       ////
///////////////////////////////////////////////////////////////// */

#ifndef HERNQUISTMODEL_HPP
#define HERNQUISTMODEL_HPP

#include "DensityModel.hpp"

//////////////////////////////////////////////////////////////////////

/** HernquistModel is a subclass of the DensityModel class and represents spherical models with a Hernquist density profile, \f[ \rho(r) = \frac{1}{2\pi}\,\frac{M_{\text{tot}}}{b^3}\left(\frac{r}{b}\right)^{-1}\left(1+ \frac{r}{b}\right)^{-3}. \f] The free parameters are the total mass \f$M_{\text{tot}}\f$ and the scale length \f$b\f$. For more information, see <a href="https://ui.adsabs.harvard.edu/abs/1990ApJ...356..359H/abstract">Hernquist (1990)</a>. */

class HernquistModel : public DensityModel
{
public:

    /** Constructor of the HernquistModel class. */
    HernquistModel(double Mtot, double b, const GaussLegendre* gl);

    /** This function returns the scale radius \f$b\f$ of the Hernquist model. */
    double scale_radius() const;
    
    /** This function returns the density \f$\rho(r)\f$ of the Hernquist model at radius \f$r\f$. */
    double density(double r) const;
    
    /** This function returns the derivative of the density \f$\rho'(r)\f$ of the Hernquist model at radius \f$r\f$. */
    double derivative_density(double r) const;
    
    /** This function returns the second derivative of the density \f$\rho''(r)\f$ of the Hernquist model at radius \f$r\f$. */
    double second_derivative_density(double r) const;
    
    /** This function returns the mass \f$M(r)\f$ of the Hernquist model at radius \f$r\f$. */
    double mass(double r) const;

    /** This function returns the total mass \f$M_{\text{tot}}\f$ of the Hernquist model. */
    double total_mass() const;

    /** This function returns the potential \f$\Psi(r)\f$ of the Hernquist model at radius \f$r\f$. */
    double potential(double r) const;
    
    /** This function returns the central potential \f$\Psi_0\f$ of the Hernquist model. */
    double central_potential() const;

private:
    
    /** The total mass \f$M_{\text{tot}}\f$. */
    double _Mtot;
    
    /** The scale radius \f$b\f$. */
    double _b;
};

//////////////////////////////////////////////////////////////////////

#endif
