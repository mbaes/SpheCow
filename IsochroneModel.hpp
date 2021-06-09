/*//////////////////////////////////////////////////////////////////
////     SpheCow -- Flexible dynamical models for galaxies      ////
////                and dark matter haloes                      ////
////     © Sterrenkundig Observatorium, Universiteit Gent       ////
///////////////////////////////////////////////////////////////// */

#ifndef ISOCHRONEMODEL_HPP
#define ISOCHRONEMODEL_HPP

#include "DensityModel.hpp"

//////////////////////////////////////////////////////////////////////

/** IsochroneModel is a subclass of the DensityModel class and represents spherical models with an isochrone density profile, \f[ \rho(r) = \frac{M_{\text{tot}}}{4\pi} \left[ \frac{3\,(b+\sqrt{r^2+b^2})\,(r^2+b^2) - r^2\,(b+3\sqrt{r^2+b^2})}{(b+\sqrt{r^2+b^2})^3\,(r^2+b^2)^{3/2}} \right]. \f] The free parameters are the total mass \f$M_{\text{tot}}\f$ and the scale length \f$b\f$. For more information, see <a href="https://ui.adsabs.harvard.edu/abs/1959AnAp...22..126H/abstract">Hénon (1959)</a>. */

class IsochroneModel : public DensityModel
{
public:
    
    /** Constructor of the IsochroneModel class. */
    IsochroneModel(double Mtot, double b, const GaussLegendre* gl);
    
    /** This function returns the scale radius \f$b\f$ of the isochrone model. */
    double scale_radius() const;
    
    /** This function returns the density \f$\rho(r)\f$ of the isochrone model at radius \f$r\f$. */
    double density(double r) const;
    
    /** This function returns the derivative of the density \f$\rho'(r)\f$ of the isochrone model at radius \f$r\f$. */
    double derivative_density(double r) const;
    
    /** This function returns the second derivative of the density \f$\rho''(r)\f$ of the isochrone model at radius \f$r\f$. */
    double second_derivative_density(double r) const;
    
    /** This function returns the mass \f$M(r)\f$ of the isochrone model at radius \f$r\f$. */
    double mass(double r) const;
    
    /** This function returns the total mass \f$M_{\text{tot}}\f$ of the isochrone model. */
    double total_mass() const;

    /** This function returns the potential \f$\Psi(r)\f$ of the isochrone model at radius \f$r\f$. */
    double potential(double r) const;

    /** This function returns the central potential \f$\Psi_0\f$ of the isochrone model. */
    double central_potential() const;

private:
    
    /** The total mass \f$M_{\text{tot}}\f$. */
    double _Mtot;
    
    /** The scale radius \f$b\f$. */
    double _b;
};

//////////////////////////////////////////////////////////////////////

#endif
