/*//////////////////////////////////////////////////////////////////
////     SpheCow -- Flexible dynamical models for galaxies      ////
////                and dark matter haloes                      ////
////     © Sterrenkundig Observatorium, Universiteit Gent       ////
///////////////////////////////////////////////////////////////// */

#ifndef SIGMOIDDENSITYMODEL_HPP
#define SIGMOIDDENSITYMODEL_HPP

#include "DensityModel.hpp"

//////////////////////////////////////////////////////////////////////

/** SigmoidDensityModel is a subclass of the DensityModel class and represents spherical models with algebraic sigmoid function as density slope, \f[ \rho(r) \propto \left(\frac{r}{r_{\text{b}}}\right)^{-\frac{\beta+\gamma}{2}} \exp\left[-\frac{\beta-\gamma}{2\alpha} \sqrt{1+\alpha^2\ln^2\left(\frac{r}{r_{\text{b}}}\right)}\right]. \f] The free parameters are the total mass \f$M_{\text{tot}}\f$, the break radius \f$r_{\text{b}}\f$, the smoothness parameter \f$\alpha\f$, the outer density slope \f$\beta\f$, and the inner density slope \f$\gamma\f$. */

class SigmoidDensityModel : public DensityModel
{
public:

    /** Constructor of the SigmoidDensityModel class. */
    SigmoidDensityModel(double Mtot, double rb, double alpha, double beta, double gamma, const GaussLegendre* gl);
    
    /** This function returns the break radius \f$r_{\text{b}}\f$ of the sigmoid density model. */
    double scale_radius() const;

    /** This function returns the density \f$\rho(r)\f$ of the sigmoid density model at radius \f$r\f$. */
    double density(double r) const;

    /** This function returns the derivative of the density \f$\rho'(r)\f$ of the sigmoid density model at radius \f$r\f$. */
    double derivative_density(double r) const;

    /** This function returns the second derivative of the density \f$\rho''(r)\f$ of the sigmoid density model at radius \f$r\f$. */
    double second_derivative_density(double r) const;

private:

    /** The total mass \f$M_{\text{tot}}\f$. */
    double _Mtot;

    /** The break radius \f$r_{\text{b}}\f$. */
    double _rb;

    /** The smoothness parameter \f$\alpha\f$. */
    double _alpha;

    /** The outer density slope \f$\beta\f$. */
    double _beta;

    /** The inner density slope \f$\gamma\f$. */
    double _gamma;
    
    /** A density pre-factor. */
    double _rhoc;
};

//////////////////////////////////////////////////////////////////////

#endif
