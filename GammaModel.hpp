/*//////////////////////////////////////////////////////////////////
////     SpheCow -- Flexible dynamical models for galaxies      ////
////                and dark matter haloes                      ////
////     © Sterrenkundig Observatorium, Universiteit Gent       ////
///////////////////////////////////////////////////////////////// */

#ifndef GAMMAMODEL_HPP
#define GAMMAMODEL_HPP

#include "DensityModel.hpp"

//////////////////////////////////////////////////////////////////////

/** GammaModel is a subclass of the DensityModel class and represents spherical models with a Dehnen or \f$\gamma\f$-density profile, \f[ \rho(r) = \frac{3-\gamma}{4\pi}\,\frac{M_{\text{tot}}}{b^3}\left(\frac{r}{b}\right)^{-\gamma}\left(1+ \frac{r}{b}\right)^{\gamma-4}. \f] The free parameters are the total mass \f$M_{\text{tot}}\f$, the scale length \f$b\f$ and the central density slope \f$\gamma\f$. For more information, see <a href="https://ui.adsabs.harvard.edu/abs/1993MNRAS.265..250D/abstract">Dehnen (1993)</a> and <a href="https://ui.adsabs.harvard.edu/abs/1994AJ....107..634T/abstract">Tremaine et al. (1994)</a>. */

class GammaModel : public DensityModel
{
public:

    /** Constructor of the GammaModel class. */
    GammaModel(double Mtot, double b, double gamma, const GaussLegendre* gl);
    
    /** This function returns the scale radius \f$b\f$ of the \f$\gamma\f$-model. */
    double scale_radius() const;

    /** This function returns the density \f$\rho(r)\f$ of the \f$\gamma\f$-model at radius \f$r\f$. */
    double density(double r) const;
    
    /** This function returns the derivative of the density \f$\rho'(r)\f$ of the \f$\gamma\f$-model at radius \f$r\f$. */
    double derivative_density(double r) const;

    /** This function returns the second derivative of the density \f$\rho''(r)\f$ of the \f$\gamma\f$-model at radius \f$r\f$. */
    double second_derivative_density(double r) const;

    /** This function returns the mass \f$M(r)\f$ of the \f$\gamma\f$-model at radius \f$r\f$. */
    double mass(double r) const;
    
    /** This function returns the potential \f$\Psi(r)\f$ of the \f$\gamma\f$-model at radius \f$r\f$. */
    double potential(double r) const;
    
private:
    
    /** The total mass \f$M_{\text{tot}}\f$. */
    double _Mtot;

    /** The scale radius \f$b\f$. */
    double _b;

    /** The central logarithmic density slope \f$\gamma\f$. */
    double _gamma;
    
    /** The density at the scale radius. */
    double _rhob;
};

//////////////////////////////////////////////////////////////////////

#endif
