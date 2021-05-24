/*//////////////////////////////////////////////////////////////////
////     SpheCow -- Flexible dynamical models for galaxies      ////
////                and dark matter haloes                      ////
////     © Sterrenkundig Observatorium, Universiteit Gent       ////
///////////////////////////////////////////////////////////////// */

#ifndef SERSICMODEL_HPP
#define SERSICMODEL_HPP

#include "SurfaceDensityModel.hpp"

//////////////////////////////////////////////////////////////////////

/** SersicModel is a subclass of the SurfaceDensityModel class and represents spherical models with a Sérsic surface density profile, \f[ \Sigma(R) =  \frac{b^{2m}}{2\pi\,m\,\Gamma(2m)}\, \frac{M_{\text{tot}}}{R_{\text{eff}}^2} \exp\left[-b\left(\frac{R}{R_{\text{eff}}}\right)^{1/m}\right]. \f] The free parameters are the total mass \f$M_{\text{tot}}\f$, the effective radius \f$R_{\text{eff}}\f$, and the Sérsic index \f$m\f$. The parameter \f$b\f$ is not a free parameter, but a numerical constant that depends on the Sérsic index. For more information, see <a href="https://ui.adsabs.harvard.edu/abs/1991A%26A...249...99C/abstract">Ciotti (1991)</a>  and <a href="https://ui.adsabs.harvard.edu/abs/2019A%26A...626A.110B/abstract">Baes & Ciotti (2020)</a>. */

class SersicModel : public SurfaceDensityModel
{
public:
    
    /** Constructor of the SersicModel class. */
    SersicModel(double Mtot, double Reff, double m, const GaussLegendre* gl);
    
    /** This function returns the effective radius \f$R_{\text{eff}}\f$ of the Sérsic model. */
    double scale_radius() const;
    
    /** This function returns the surface density \f$\Sigma(R)\f$ of the Sérsic model at projected radius \f$R\f$. */
    double surface_density(double R) const;

    /** This function returns the derivative of the surface density \f$\Sigma'(R)\f$ of the Sérsic model at projected radius \f$R\f$. */
    double derivative_surface_density(double R) const;

    /** This function returns the second derivative of the surface density \f$\Sigma''(R)\f$ of the Sérsic model at projected radius \f$R\f$. */
    double second_derivative_surface_density(double R) const;
    
    /** This function returns the third derivative of the surface density \f$\Sigma'''(R)\f$ of the Sérsic model at projected radius \f$R\f$. */
    double third_derivative_surface_density(double R) const;
    
private:
    
    /** The total mass \f$M_{\text{tot}}\f$. */
    double _Mtot;
    
    /** The effective radius \f$R_{\text{eff}}\f$. */
    double _Reff;
    
    /** The Sérsic index \f$m\f$. */
    double _m;
    
    /** The numerical constant \f$b\f$. */
    double _b;

    /** The central surface brightness. */
    double _Sigma0;
};


//////////////////////////////////////////////////////////////////////

#endif
