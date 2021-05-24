/*//////////////////////////////////////////////////////////////////
////     SpheCow -- Flexible dynamical models for galaxies      ////
////                and dark matter haloes                      ////
////     © Sterrenkundig Observatorium, Universiteit Gent       ////
///////////////////////////////////////////////////////////////// */

#ifndef PERFECTSPHEREMODEL_HPP
#define PERFECTSPHEREMODEL_HPP

#include "DensityModel.hpp"

//////////////////////////////////////////////////////////////////////

/** PerfectSphereModel is a subclass of the DensityModel class and represents spherical models with a perfect sphere density profile, \f[ \rho(r) = \frac{1}{\pi^2}\,\frac{M_{\text{tot}}}{c^3} \left(1+ \frac{r^2}{c^2}\right)^{-2}. \f] The free parameters are the total mass \f$M_{\text{tot}}\f$ and the scale length \f$c\f$. For more information, see <a href="https://ui.adsabs.harvard.edu/abs/1985MNRAS.216..273D/abstract">de Zeeuw (1985)</a>. */

class PerfectSphereModel : public DensityModel
{
public:

    /** Constructor of the PerfectSphereModel class. */
    PerfectSphereModel(double Mtot, double c, const GaussLegendre* gl);
    
    /** This function returns the scale radius \f$c\f$ of the perfect sphere model. */
    double scale_radius() const;

    /** This function returns the density \f$\rho(r)\f$ of the perfect sphere model at radius \f$r\f$. */
    double density(double r) const;

    /** This function returns the derivative of the density \f$\rho'(r)\f$ of the perfect sphere model at radius \f$r\f$. */
    double derivative_density(double r) const;

    /** This function returns the second derivative of the density \f$\rho''(r)\f$ of the perfect sphere model at radius \f$r\f$. */
    double second_derivative_density(double r) const;

    /** This function returns the mass \f$M(r)\f$ of the perfect sphere model at radius \f$r\f$. */
    double mass(double r) const;

    /** This function returns the potential \f$\Psi(r)\f$ of the perfect sphere model at radius \f$r\f$. */
    double potential(double r) const;

private:

    /** The total mass \f$M_{\text{tot}}\f$. */
    double _Mtot;

    /** The scale radius \f$c\f$. */
    double _c;
};

//////////////////////////////////////////////////////////////////////

#endif
