/*//////////////////////////////////////////////////////////////////
////     SpheCow -- Flexible dynamical models for galaxies      ////
////                and dark matter haloes                      ////
////     © Sterrenkundig Observatorium, Universiteit Gent       ////
///////////////////////////////////////////////////////////////// */

#ifndef BPLMODEL_HPP
#define BPLMODEL_HPP

#include "DensityModel.hpp"

//////////////////////////////////////////////////////////////////////

/** BPLModel is a subclass of the DensityModel class and represents spherical models with a broken power-law density profile, \f[ \rho(r) = \frac{(\beta-3)\,(3-\gamma)}{4\pi\,(\beta-\gamma)}\, \frac{M_{\text{tot}}}{r_{\text{b}}^3} \times \begin{cases} \displaystyle \; \left(\frac{r}{r_{\text{b}}}\right)^{-\gamma} & r\leq r_{\text{b}}, \\ \displaystyle \; \left(\frac{r}{r_{\text{b}}}\right)^{-\beta} & r\geq r_{\text{b}}. \end{cases} \f] The free parameters are the total mass \f$M_{\text{tot}}\f$, the break radius \f$r_{\text{b}}\f$, the outer density slope \f$\beta\f$, and the inner density slope \f$\gamma\f$. For more information, see <a href="https://ui.adsabs.harvard.edu/abs/2021MNRAS.503.2955B/abstract">Baes & Camps (2021)</a>. */

class BPLModel : public DensityModel
{
public:
    
    /** Constructor of the BPLModel class. */
    BPLModel(double Mtot, double rb, double beta, double gamma, const GaussLegendre* gl);
    
    /** This function returns the break radius \f$r_{\text{b}}\f$ of the BPL model. */
    double scale_radius() const;
    
    /** This function returns the density \f$\rho(r)\f$ of the BPL model at radius \f$r\f$. */
    double density(double r) const;
    
    /** This function returns the derivative of the density \f$\rho'(r)\f$ of the BPL model at radius \f$r\f$. */
    double derivative_density(double r) const;
    
    /** This function returns the second derivative of the density \f$\rho''(r)\f$ of the BPL model at radius \f$r\f$. */
    double second_derivative_density(double r) const;
    
    /** This function returns the mass \f$M(r)\f$ of the BPL model at radius \f$r\f$. */
    double mass(double r) const;
    
    /** This function returns the potential \f$\Psi(r)\f$ of the BPL model at radius \f$r\f$. */
    double potential(double r) const;
    
    /** This function returns the isotropic distribution function \f$f_{\text{iso}}({\cal{E}})\f$ of the BPL model at radius \f$r=r(\cal{E})\f$. Since the second derivative of the density is discontinuous, the general Eddington formula needs to be complimented with an additional term. */
    double isotropic_distribution_function(double r) const;

    /** This function returns the Osipkov-Merritt distribution function \f$f_{\text{om}}(Q)\f$ of the BPL model at radius \f$r=r(Q)\f$, for an anisotropy radius \f$r_{\text{a}}\f$. Since the second derivative of the density is discontinuous, the general Osipkov-Merritt formula needs to be complimented with an additional term. */
    double osipkov_merritt_distribution_function(double r, double ra) const;
    
private:
    
    /** The total mass \f$M_{\text{tot}}\f$. */
    double _Mtot;
    
    /** The break radius \f$r_{\text{b}}\f$. */
    double _rb;
    
    /** The outer density slope \f$\beta\f$. */
    double _beta;

    /** The inner density slope \f$\beta\f$. */
    double _gamma;
    
    /** A density pre-factor. */
    double _rhoff;
};

//////////////////////////////////////////////////////////////////////

#endif
