/*//////////////////////////////////////////////////////////////////
////     SpheCow -- Flexible dynamical models for galaxies      ////
////                and dark matter haloes                      ////
////     © Sterrenkundig Observatorium, Universiteit Gent       ////
///////////////////////////////////////////////////////////////// */

#ifndef GAUSSLEGENDRE_HPP
#define GAUSSLEGENDRE_HPP

#include "Basics.hpp"
#include <functional>

//////////////////////////////////////////////////////////////////////

/** GaussLegendre is the class that represents Gauss-Legendre integrators. A Gauss-Legendre integrator approximates a Riemann integral by a finite sum \f[ \int_0^1 f(x)\,{\text{d}}x \approx \sum_{i=1}^N w_i\,f(x_i). \f] For a given value of \f$N\f$, the nodes \f$x_i\f$ and the weights \f$w_i\f$ are fixed, and correspond to the unique choice that allows the quadrature rule to integrate all polynomials of degree \f$2N+1\f$ exactly. For more background information on Gauss-Legendre quadrature, see the <a href="https://en.wikipedia.org/wiki/Gauss–Legendre_quadrature">Wikipedia page</a>. Note that the standard convention for Gauss-Legendre integration is the interval \f$[-1,1]\f$, but we adopt the interval \f$[0,1]\f$ because it serves our purposes better. */

class GaussLegendre
{
public:
    
    /** Constructor for the GaussLegendre class. The constructor reads in a number \f$n\f$ and constructs a Gauss-Legendre integrator with \f$N\f$ nodes, where \f$N\f$ is the smallest power of 2 that is equal to or larger than \f$n\f$. The constructor just reads the appropriate list of nodes and weigths from a precalculated table. The minimum number of nodes is 8, the maximum number is 512. */
    GaussLegendre(int num);
    
    /** This function returns an estimate of the integral of the function \f$X(u)\f$ over the interval \f$[0,+\infty[\f$. The integral is split at the break radius \f$r_{\text{b}}\f$, and converted to \f[ \int_0^\infty X(u)\, {\text{d}}u = r_{\text{b}} \int_0^{\pi/2} X(r_{\text{b}}\sin\theta) \cos\theta\, {\text{d}}\theta + r_{\text{b}} \int_0^{\pi/2} X(r_{\text{b}}\csc\theta) \cos\theta\, \csc^2\theta\, {\text{d}}\theta. \f] The resulting integrals are estimated as Gauss-Legendre quadratures. */
    double integrate_0_infty(std::function<double(double)> X, double rb) const;
    
    /** This function returns an estimate of the integral of the function \f$X(u)\f$ over the interval \f$[0,r]\f$, with \f$r>0\f$ an arbitrary number. If \f$r\leq r_{\text{b}}\f$, the integral is converted to \f[ \int_0^r X(u)\, {\text{d}}u = r \int_0^{\pi/2} X(r \sin\theta) \cos\theta\, {\text{d}}\theta. \f] If \f$r>r_{\text{b}}\f$, the integral is split at the break radius \f$r_{\text{b}}\f$, and converted to \f[ \int_0^r X(u)\, {\text{d}}u = r_{\text{b}} \int_0^{\pi/2} X(r_{\text{b}}\sin\theta) \cos\theta\, {\text{d}}\theta + r \int_{\arcsin(r_{\text{b}}/r)}^{\pi/2} X(r\sin\theta) \cos\theta\, {\text{d}}\theta. \f] The resulting integrals are estimated as Gauss-Legendre quadratures. */
    double integrate_0_r(std::function<double(double)> X, double r, double rb) const;
    
    /** This function returns an estimate of the integral of the function \f$X(u)\f$ over the interval \f$[r,+\infty[\f$, with \f$r>0\f$ an arbitrary number. If \f$r<r_{\text{b}}\f$, the integral is split at the break radius \f$r_{\text{b}}\f$, and converted to \f[ \int_r^\infty X(u)\, {\text{d}}u = r \int_{\arcsin(r/r_{\text{b}})}^{\pi/2} X(r \csc\theta) \cos\theta \csc^2\theta\, {\text{d}}\theta + r_{\text{b}} \int_0^{\pi/2} X(r_{\text{b}}\csc\theta) \cos\theta \csc^2\theta\, {\text{d}}\theta. \f] If \f$r\geq r_{\text{b}}\f$, the integral is converted to \f[ \int_0^r X(u)\, {\text{d}}u = r \int_0^{\pi/2} X(r \csc\theta) \cos\theta \csc^2\theta\, {\text{d}}\theta. \f] The resulting integrals are estimated as Gauss-Legendre quadratures. */
    double integrate_r_infty(std::function<double(double)> X, double r, double rb) const;
    
private:
    
    /** The number of nodes \f$N\f$. */
    int _num;

    /** A vector with the \f$N\f$ nodes \f$x_i\f$. */
    std::vector<double> _xv;

    /** A vector with the \f$N\f$ weights \f$w_i\f$. */
    std::vector<double> _wv;
};

//////////////////////////////////////////////////////////////////////

#endif
