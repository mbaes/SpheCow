/*//////////////////////////////////////////////////////////////////
////     SpheCow -- Flexible dynamical models for galaxies      ////
////                and dark matter haloes                      ////
////     © Sterrenkundig Observatorium, Universiteit Gent       ////
///////////////////////////////////////////////////////////////// */

#ifndef PLUMMERMODEL_HPP
#define PLUMMERMODEL_HPP

#include "DensityModel.hpp"

//////////////////////////////////////////////////////////////////////

/** PlummerModel is a subclass of the DensityModel class and represents spherical models with a Plummer density profile, \f[ \rho(r) = \frac{3}{4\pi}\,\frac{M_{\text{tot}}}{c^3} \left(1+ \frac{r^2}{c^2}\right)^{-5/2}. \f] The free parameters are the total mass \f$M_{\text{tot}}\f$ and the scale length \f$c\f$. For more information, see <a href="https://ui.adsabs.harvard.edu/abs/1987MNRAS.224...13D/abstract">Dejonghe (1987)</a>. */

class PlummerModel : public DensityModel
{
public:

    /** Constructor of the PlummerModel class. */
    PlummerModel(double Mtot, double c, const GaussLegendre* gl);

    /** This function returns the scale radius \f$c\f$ of the Plummer model. */
    double scale_radius() const;

    /** This function returns the density \f$\rho(r)\f$ of the Plummer model at radius \f$r\f$. */
    double density(double r) const;

    /** This function returns the derivative of the density \f$\rho'(r)\f$ of the Plummer model at radius \f$r\f$. */
    double derivative_density(double r) const;

    /** This function returns the second derivative of the density \f$\rho''(r)\f$ of the Plummer model at radius \f$r\f$. */
    double second_derivative_density(double r) const;

    /** This function returns the mass \f$M(r)\f$ of the Plummer model at radius \f$r\f$. */
    double mass(double r) const;

    /** This function returns the potential \f$\Psi(r)\f$ of the Plummer model at radius \f$r\f$. */
    double potential(double r) const;

private:
    
    /** The total mass \f$M_{\text{tot}}\f$. */
    double _Mtot;

    /** The scale radius \f$c\f$. */
    double _c;
};

//////////////////////////////////////////////////////////////////////

#endif
