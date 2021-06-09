/*//////////////////////////////////////////////////////////////////
////     SpheCow -- Flexible dynamical models for galaxies      ////
////                and dark matter haloes                      ////
////     © Sterrenkundig Observatorium, Universiteit Gent       ////
///////////////////////////////////////////////////////////////// */

#ifndef DEVAUCOULEURSMODEL_HPP
#define DEVAUCOULEURSMODEL_HPP

#include "SurfaceDensityModel.hpp"

//////////////////////////////////////////////////////////////////////

/** DeVaucouleursModel is a subclass of the SurfaceDensityModel class and represents spherical models with a de Vaucouleurs surface density profile, \f[ \Sigma(R) =  \frac{b_4^{8}}{40320\pi}\,\frac{M_{\text{tot}}}{R_{\text{eff}}^2}\exp\left[-b_4\left(\frac{R}{R_{\text{eff}}}\right)^{1/4}\right], \f] with \f$b_4 = 7.669249443\f$ a numerical constant. The free parameters are the total mass \f$M_{\text{tot}}\f$ and the effective radius \f$R_{\text{eff}}\f$. For more information, see <a href="https://ui.adsabs.harvard.edu/abs/1991A%26A...249...99C/abstract">Ciotti (1991)</a>  and <a href="https://ui.adsabs.harvard.edu/abs/2019A%26A...626A.110B/abstract">Baes & Ciotti (2020)</a>. */

class DeVaucouleursModel : public SurfaceDensityModel
{
public:
    
    /** Constructor of the DeVaucouleursModel class. */
    DeVaucouleursModel(double Mtot, double Reff, const GaussLegendre* gl);
    
    /** This function returns the effective radius \f$R_{\text{eff}}\f$ of the de Vaucouleurs model. */
    double scale_radius() const;
    
    /** This function returns the surface density \f$\Sigma(R)\f$ of the de Vaucouleurs model at projected radius \f$R\f$. */
    double surface_density(double R) const;

    /** This function returns the derivative of the surface density \f$\Sigma'(R)\f$ of the de Vaucouleurs model at projected radius \f$R\f$. */
    double derivative_surface_density(double R) const;
    
    /** This function returns the second derivative of the surface density \f$\Sigma''(R)\f$ of the de Vaucouleurs model at projected radius \f$R\f$. */
    double second_derivative_surface_density(double R) const;
    
    /** This function returns the third derivative of the surface density \f$\Sigma'''(R)\f$ of the de Vaucouleurs model at projected radius \f$R\f$. */
    double third_derivative_surface_density(double R) const;
    
    /** This function returns the total mass \f$M_{\text{tot}}\f$ of the de Vaucouleurs model. */
    double total_mass() const;

    /** This function returns the central potential \f$\Psi_0\f$ of the de Vaucouleurs model. */
    double central_potential() const;

private:
    
    /** The total mass \f$M_{\text{tot}}\f$. */
    double _Mtot;
    
    /** The effective radius \f$R_{\text{eff}}\f$. */
    double _Reff;
    
    /** The numerical constant \f$b_4\f$. */
    double _b;
    
    /** The pre-factor for the surface density. */
    double _Sigmaff;
};


//////////////////////////////////////////////////////////////////////

#endif
