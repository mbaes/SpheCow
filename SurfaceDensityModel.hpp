/*//////////////////////////////////////////////////////////////////
////     SpheCow -- Flexible dynamical models for galaxies      ////
////                and dark matter haloes                      ////
////     © Sterrenkundig Observatorium, Universiteit Gent       ////
///////////////////////////////////////////////////////////////// */

#ifndef SURFACEDENSITYMODEL_HPP
#define SURFACEDENSITYMODEL_HPP

#include "Model.hpp"

//////////////////////////////////////////////////////////////////////

/** SurfaceDensityModel is an abstract subclass of the Model class and represents the base class for spherical models defined through their surface density profile. */

class SurfaceDensityModel : public Model
{
public:

    /** Virtual destructor of the SurfaceDensityModel class. */
    virtual ~SurfaceDensityModel() {};
    
    /** This pure virtual function returns the second derivative of the surface density \f$\Sigma''(R)\f$ at projected radius \f$R\f$. */
    virtual double second_derivative_surface_density(double R) const = 0;
    
    /** This pure virtual function returns the third derivative of the surface density \f$\Sigma'''(R)\f$ at projected radius \f$R\f$. */
    virtual double third_derivative_surface_density(double R) const = 0;
    
    /** This function returns the density \f$\rho(r)\f$ at radius \f$r\f$. It is calculated as \f[ \rho(r) = -\frac{1}{\pi} \int_r^\infty \frac{\Sigma'(u)\,{\text{d}} u}{\sqrt{u^2-r^2}}. \f] The integration is performed using Gauss-Legendre quadrature. */
    double density(double r) const;
    
    /** This function returns the derivative of the density \f$\rho'(r)\f$ at radius \f$r\f$.  It is calculated as \f[  \rho'(r) = -\frac{1}{\pi} \int_r^\infty \frac{\Sigma''(u)\,u\,{\text{d}} u}{r\sqrt{u^2-r^2}}. \f] The integration is performed using Gauss-Legendre quadrature.  */
    double derivative_density(double r) const;
    
    /** This function returns the second derivative of the density \f$\rho''(r)\f$ at radius \f$r\f$.  It is calculated as \f[  \rho''(r) = -\frac{1}{\pi} \int_r^\infty \frac{\Sigma'''(u)\,u^2\,{\text{d}}u}{r^2\sqrt{u^2-r^2}}. \f] The integration is performed using Gauss-Legendre quadrature.  */
    double second_derivative_density(double r) const;
    
    /** This function returns the mass \f$M(r)\f$ at radius \f$r\f$. It is calculated as \f[ M(r) = -\pi \left[ \int_0^r \Sigma'(u)\,u^2\, {\text{d}} u + \int_r^\infty \Sigma'(u)\,w_-(u,r)\, {\text{d}} u \right],\f] with \f[ w_-(u,r) = \frac{2}{\pi}\left[u^2\arctan\left(\frac{r}{\sqrt{u^2-r^2}}\right)-r\sqrt{u^2-r^2}\right]. \f] The integration is performed using Gauss-Legendre quadrature. */
    double mass(double r) const;
    
    /** This function returns the total mass \f$M_{\text{tot}}\f$. It is calculated as \f[ M_{\text{tot}} = 2\pi \int_0^\infty \Sigma(u)\, u\, {\text{d}} u. \f] The integration is performed using Gauss-Legendre quadrature.  */
    double total_mass() const;
    
    /** This function returns the potential \f$\Psi(r)\f$ at radius \f$r\f$. It is calculated as \f[ \Psi(r) = -\frac{\pi\,G}{r} \left[ \int_0^r \Sigma'(u)\,u^2\, {\text{d}} u + \int_r^\infty \Sigma'(u)\,w_+(u,r)\,{\text{d}} u \right], \f] with \f[ w_+(u,r) = \frac{2}{\pi}\left[u^2\arctan\left(\frac{r}{\sqrt{u^2-r^2}}\right)+r\sqrt{u^2-r^2}\right]. \f] The integration is performed using Gauss-Legendre quadrature. */
    double potential(double r) const;
};

//////////////////////////////////////////////////////////////////////

#endif
