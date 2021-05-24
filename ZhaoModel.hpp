/*//////////////////////////////////////////////////////////////////
////     SpheCow -- Flexible dynamical models for galaxies      ////
////                and dark matter haloes                      ////
////     © Sterrenkundig Observatorium, Universiteit Gent       ////
///////////////////////////////////////////////////////////////// */

#ifndef ZHAOMODEL_HPP
#define ZHAOMODEL_HPP

#include "DensityModel.hpp"

//////////////////////////////////////////////////////////////////////

/** ZhaoModel is a subclass of the DensityModel class and represents spherical models with a Zhao density profile, \f[ \rho(r) =  \frac{\alpha}{4\pi}\, \frac{\Gamma\left(\frac{\beta-\gamma}{\alpha}\right)}{\Gamma\left(\frac{\beta-3}{\alpha}\right)\, \Gamma\left(\frac{3-\gamma}{\alpha}\right)}\, \frac{M_{\text{tot}}}{r_{\text{b}}^3}\, \left(\frac{r}{r_{\text{b}}}\right)^{-\gamma} \left[ 1+ \left(\frac{r}{r_{\text{b}}}\right)^\alpha\right]^{\frac{\gamma-\beta}{\alpha}}. \f] The free parameters are the total mass \f$M_{\text{tot}}\f$, the break radius \f$r_{\text{b}}\f$, the smoothness parameter \f$\alpha\f$, the outer density slope \f$\beta\f$, and the inner density slope \f$\gamma\f$. For more information, see <a href="https://ui.adsabs.harvard.edu/abs/1996MNRAS.278..488Z/abstract">Zhao (1996)</a>. Note that we use a different convention for \f$\alpha\f$ than <a href="https://ui.adsabs.harvard.edu/abs/1996MNRAS.278..488Z/abstract">Zhao (1996)</a>: in our case, larger values of \f$\alpha\f$ correspond to sharper breaks between the inner and outer density profiles. */

class ZhaoModel : public DensityModel
{
public:

    /** Constructor of the ZhaoModel class. */
    ZhaoModel(double Mtot, double rb, double alpha, double beta, double gamma, const GaussLegendre* gl);

    /** This function returns the break radius \f$r_{\text{b}}\f$ of the Zhao model. */
    double scale_radius() const;
    
    /** This function returns the density \f$\rho(r)\f$ of the Zhao model at radius \f$r\f$. */
    double density(double r) const;
    
    /** This function returns the derivative of the density \f$\rho'(r)\f$ of the Zhao model at radius \f$r\f$. */
    double derivative_density(double r) const;
    
    /** This function returns the second derivative of the density \f$\rho''(r)\f$ of the Zhao model at radius \f$r\f$. */
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

    /** The inner density slope \f$\beta\f$. */
    double _gamma;
    
    /** A density pre-factor. */
    double _rhoff;
};

//////////////////////////////////////////////////////////////////////

#endif
