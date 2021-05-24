/*//////////////////////////////////////////////////////////////////
////     SpheCow -- Flexible dynamical models for galaxies      ////
////                and dark matter haloes                      ////
////     © Sterrenkundig Observatorium, Universiteit Gent       ////
///////////////////////////////////////////////////////////////// */

#ifndef HYPERVIRIALMODEL_HPP
#define HYPERVIRIALMODEL_HPP

#include "DensityModel.hpp"

//////////////////////////////////////////////////////////////////////

/** HypervirialModel is a subclass of the DensityModel class and represents spherical models with a hypervirial density profile, \f[ \rho(r) = \frac{p+1}{4\pi}\, \frac{M_{\text{tot}}}{r_{\text{s}}^3} \left(\frac{r}{r_{\text{s}}}\right)^{p-2}\left[1+ \left(\frac{r}{r_{\text{s}}}\right)^p\right]^{-2-1/p}. \f] The free parameters are the total mass \f$M_{\text{tot}}\f$, the scale length \f$r_{\text{s}}\f$, and the hypervirial index \f$p\f$. For more information, see <a href="https://ui.adsabs.harvard.edu/abs/2005MNRAS.360..492E/abstract">Evans & An (2005)</a>. */

class HypervirialModel : public DensityModel
{
public:
    /** Constructor of the HypervirialModel class. */
    HypervirialModel(double Mtot, double rs, double p, const GaussLegendre* gl);
    
    /** This function returns the scale radius \f$r_{\text{s}}\f$ of the hypervirial model. */
    double scale_radius() const;
    
    /** This function returns the density \f$\rho(r)\f$ of the hypervirial model at radius \f$r\f$. */
    double density(double r) const;
    
    /** This function returns the derivative of the density \f$\rho'(r)\f$ of the hypervirial model at radius \f$r\f$. */
    double derivative_density(double r) const;
    
    /** This function returns the second derivative of the density \f$\rho''(r)\f$ of the hypervirial model at radius \f$r\f$. */
    double second_derivative_density(double r) const;
    
    /** This function returns the mass \f$M(r)\f$ of the hypervirial model at radius \f$r\f$. */
    double mass(double r) const;
    
    /** This function returns the potential \f$\Psi(r)\f$ of the hypervirial model at radius \f$r\f$. */
    double potential(double r) const;

private:
    
    /** The total mass \f$M_{\text{tot}}\f$. */
    double _Mtot;
    
    /** The scale radius \f$r_{\text{s}}\f$. */
    double _rs;
    
    /** The hypervirial index \f$p\f$.*/
    double _p;
};

//////////////////////////////////////////////////////////////////////

#endif
