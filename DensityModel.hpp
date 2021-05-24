/*//////////////////////////////////////////////////////////////////
////     SpheCow -- Flexible dynamical models for galaxies      ////
////                and dark matter haloes                      ////
////     © Sterrenkundig Observatorium, Universiteit Gent       ////
///////////////////////////////////////////////////////////////// */

#ifndef DENSITYMODEL_HPP
#define DENSITYMODEL_HPP

#include "Model.hpp"

//////////////////////////////////////////////////////////////////////

/** DensityModel is an abstract subclass of the Model class and represents the base class for spherical models defined through their density profile. */

class DensityModel : public Model
{
public:

    /** Virtual destructor of the DensityModel class. */
    virtual ~DensityModel() {};
    
    /** This function returns the mass \f$M(r)\f$ at radius \f$r\f$. It is calculated as \f[ M(r) = 4\pi \int_0^r \rho(u)\, u^2\, {\text{d}}u. \f] The integration is performed using Gauss-Legendre quadrature. This function is a virtual function that can be reimplemented by derived classes. */
    virtual double mass(double r) const;

    /** This function returns the total mass \f$M_{\text{tot}}\f$ . It is calculated as \f[ M_{\text{tot}} = 4\pi \int_0^\infty \rho(u)\, u^2\, {\text{d}} u. \f] The integration is performed using Gauss-Legendre quadrature. */
    double total_mass() const;
    
    /** This function returns the potential \f$\Psi(r)\f$ at radius \f$r\f$. It is calculated as \f[ \Psi(r) = \frac{GM(r)}{r} + 4\pi\,G\int_r^\infty \rho(u)\,u\,{\text{d}} u. \f] The integration is performed using Gauss-Legendre quadrature. This function is a virtual function that can be reimplemented by derived classes. */
    virtual double potential(double r) const;
    
    /** This function returns the surface density \f$\Sigma(R)\f$ at projected radius \f$R\f$. It is calculated as \f[ \Sigma(R) = 2\int_R^\infty \frac{\rho(u)\,u\,{\text{d}} u}{\sqrt{u^2-R^2}}. \f] */
    double surface_density(double R) const;
    
    /** This function returns the derivative of the surface density \f$\Sigma'(R)\f$ at projected radius \f$R\f$. It is calculated as \f[ \Sigma'(R) = 2\int_R^\infty \frac{[\rho(u)+ u\,\rho'(u)]\,u\,{\text{d}} u}{R \sqrt{u^2-R^2}}. \f] */
    double derivative_surface_density(double R) const;
};

//////////////////////////////////////////////////////////////////////

#endif
