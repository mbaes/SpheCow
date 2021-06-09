/*//////////////////////////////////////////////////////////////////
////     SpheCow -- Flexible dynamical models for galaxies      ////
////                and dark matter haloes                      ////
////     © Sterrenkundig Observatorium, Universiteit Gent       ////
///////////////////////////////////////////////////////////////// */

#ifndef EINASTOMODEL_HPP
#define EINASTOMODEL_HPP

#include "DensityModel.hpp"

//////////////////////////////////////////////////////////////////////

/** EinastoModel is a subclass of the DensityModel class and represents spherical models with an Einasto density profile, \f[ \rho(r) =  \frac{d^{3n}}{4\pi\,n\,\Gamma(3n)}\,\frac{M}{r_{\text{h}}^3}\exp\left[-d\left(\frac{r}{r_{\text{h}}}\right)^{1/n}\right]. \f] The free parameters are the total mass \f$M_{\text{tot}}\f$, the half-mass radius \f$r_{\text{h}}\f$, and the Einasto index \f$n\f$. The parameter \f$d\f$ is not a free parameter, but a numerical constant that depends on the Einasto index. For more information, see <a href="https://ui.adsabs.harvard.edu/abs/2012A%26A...540A..70R/abstract">Retana-Montenegro et al. (2012)</a>. */

class EinastoModel : public DensityModel
{
public:
    
    /** Constructor of the EinastoModel class. */
    EinastoModel(double Mtot, double rh, double n, const GaussLegendre* gl);
    
    /** This function returns the half-mass radius \f$r_{\text{h}}\f$ of the Einasto model. */
    double scale_radius() const;
    
    /** This function returns the density \f$\rho(r)\f$ of the Einasto model at radius \f$r\f$. */
    double density(double r) const;
    
    /** This function returns the derivative of the density \f$\rho'(r)\f$ of the Einasto model at radius \f$r\f$. */
    double derivative_density(double r) const;

    /** This function returns the second derivative of the density \f$\rho''(r)\f$ of the Einasto model at radius \f$r\f$. */
    double second_derivative_density(double r) const;

    /** This function returns the total mass \f$M_{\text{tot}}\f$ of the Einasto model. */
    double total_mass() const;

    /** This function returns the central potential \f$\Psi_0\f$ of the Einasto model. */
    double central_potential() const;

private:
    
    /** The total mass \f$M_{\text{tot}}\f$. */
    double _Mtot;

    /** The half-mass radius \f$r_{\text{h}}\f$. */
    double _rh;
    
    /** The Einasto index \f$n\f$. */
    double _n;

    /** The dimensionless constant \f$d\f$. */
    double _d;
    
    /** The central density. */
    double _rho0;
};

//////////////////////////////////////////////////////////////////////

#endif
