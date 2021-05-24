/*//////////////////////////////////////////////////////////////////
////     SpheCow -- Flexible dynamical models for galaxies      ////
////                and dark matter haloes                      ////
////     © Sterrenkundig Observatorium, Universiteit Gent       ////
///////////////////////////////////////////////////////////////// */

#ifndef SIGMOIDSURFACEDENSITYMODEL_HPP
#define SIGMOIDSURFACEDENSITYMODEL_HPP

#include "SurfaceDensityModel.hpp"

//////////////////////////////////////////////////////////////////////

/** SigmoidSurfaceDensityModel is a subclass of the SurfaceDensityModel class and represents spherical models with algebraic sigmoid function as surface density slope, \f[ \Sigma(R) \propto \left(\frac{R}{R_{\text{b}}}\right)^{-\frac{\beta+\gamma}{2}} \exp\left[-\frac{\beta-\gamma}{2\alpha} \sqrt{1+\alpha^2\ln^2\left(\frac{R}{R_{\text{b}}}\right)}\right]. \f] The free parameters are the total mass \f$M_{\text{tot}}\f$, the break radius \f$R_{\text{b}}\f$, the smoothness parameter \f$\alpha\f$, the outer surface density slope \f$\beta\f$, and the inner surface density slope \f$\gamma\f$. */

class SigmoidSurfaceDensityModel : public SurfaceDensityModel
{
public:

    /** Constructor of the SigmoidSurfaceDensityModel class. */
    SigmoidSurfaceDensityModel(double Mtot, double Rb, double alpha, double beta, double gamma, const GaussLegendre* gl);

    /** This function returns the break radius \f$R_{\text{b}}\f$ of the sigmoid surface density model. */
    double scale_radius() const;
    
    /** This function returns the surface density \f$\Sigma(R)\f$ of the sigmoid surface density model at projected radius \f$R\f$. */
    double surface_density(double R) const;

    /** This function returns the derivative of the surface density \f$\Sigma'(R)\f$ of the sigmoid surface density model at projected radius \f$R\f$. */
    double derivative_surface_density(double R) const;
    
    /** This function returns the second derivative of the surface density \f$\Sigma''(R)\f$ of the sigmoid surface density model at projected radius \f$R\f$. */
    double second_derivative_surface_density(double R) const;
    
    /** This function returns the third derivative of the surface density \f$\Sigma'''(R)\f$ of the sigmoid surface density model at projected radius \f$R\f$. */
    double third_derivative_surface_density(double R) const;

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
    
    /** A surface density pre-factor. */
    double _Sigmac;
};

//////////////////////////////////////////////////////////////////////

#endif
