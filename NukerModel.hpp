/*//////////////////////////////////////////////////////////////////
////     SpheCow -- Flexible dynamical models for galaxies      ////
////                and dark matter haloes                      ////
////     © Sterrenkundig Observatorium, Universiteit Gent       ////
///////////////////////////////////////////////////////////////// */

#ifndef NUKERMODEL_HPP
#define NUKERMODEL_HPP

#include "SurfaceDensityModel.hpp"

//////////////////////////////////////////////////////////////////////

/** NukerModel is a subclass of the SurfaceDensityModel class and represents spherical models with a Nuker surface density profile, \f[ \Sigma(R) =  \frac{\alpha}{2\pi}\, \frac{\Gamma\left(\frac{\beta-\gamma}{\alpha}\right)}{\Gamma\left(\frac{\beta-2}{\alpha}\right)\, \Gamma\left(\frac{2-\gamma}{\alpha}\right)}\, \frac{M_{\text{tot}}}{R_{\text{b}}^2}\, \left(\frac{R}{R_{\text{b}}}\right)^{-\gamma} \left[ 1+ \left(\frac{R}{R_{\text{b}}}\right)^\alpha\right]^{\frac{\gamma-\beta}{\alpha}}. \f] The free parameters are the total mass \f$M_{\text{tot}}\f$, the break radius \f$R_{\text{b}}\f$, the smoothness parameter \f$\alpha\f$, the outer surface density slope \f$\beta\f$, and the inner surface density slope \f$\gamma\f$. For more information, see <a href="https://ui.adsabs.harvard.edu/abs/2020A%26A...634A.109B/abstract">Baes (2020). */

class NukerModel : public SurfaceDensityModel
{
public:

    /** Constructor of the NukerModel class. */
    NukerModel(double Mtot, double Rb, double alpha, double beta, double gamma, const GaussLegendre* gl);
    
    /** This function returns the break radius \f$R_{\text{b}}\f$ of the Nuker model. */
    double scale_radius() const;
    
    /** This function returns the surface density \f$\Sigma(R)\f$ of the Nuker model at projected radius \f$R\f$. */
    double surface_density(double R) const;

    /** This function returns the derivative of the surface density \f$\Sigma'(R)\f$ of the Nuker model at projected radius \f$R\f$. */
    double derivative_surface_density(double R) const;
    
    /** This function returns the second derivative of the surface density \f$\Sigma''(R)\f$ of the Nuker model at projected radius \f$R\f$. */
    double second_derivative_surface_density(double R) const;
    
    /** This function returns the third derivative of the surface density \f$\Sigma'''(R)\f$ of the Nuker model at projected radius \f$R\f$. */
    double third_derivative_surface_density(double R) const;

    /** This function returns the total mass \f$M_{\text{tot}}\f$ of the Nuker model. */
    double total_mass() const;

    /** This function returns the central potential \f$\Psi_0\f$ of the Nuker model. */
    double central_potential() const;

private:
    
    /** The total mass \f$M_{\text{tot}}\f$. */
    double _Mtot;

    /** The break radius \f$R_{\text{b}}\f$. */
    double _Rb;

    /** The smoothness parameter \f$\alpha\f$. */
    double _alpha;

    /** The outer surface density slope \f$\beta\f$. */
    double _beta;

    /** The inner surface density slope \f$\gamma\f$. */
    double _gamma;
    
    /** The surface density at the break radius. */
    double _Sigmab;
};

//////////////////////////////////////////////////////////////////////

#endif
