/*//////////////////////////////////////////////////////////////////
////     SpheCow -- Flexible dynamical models for galaxies      ////
////                and dark matter haloes                      ////
////     © Sterrenkundig Observatorium, Universiteit Gent       ////
///////////////////////////////////////////////////////////////// */

#ifndef SPHECOW_HPP
#define SPHECOW_HPP

#include "Basics.hpp"

class Model;

//////////////////////////////////////////////////////////////////////

/** This routine is the main workhorse of the SpheCow code. It calculates the most important photometric and dynamical properties for a given model, both for an isotropic orbital structure and an Osipkov-Merritt orbital structure with anisotropy radius \f$r_{\text{a}}\f$. The routine reads in a vector with radii, calculates the entire set of properties at each of these radii, and writes the results to a file.  */
void run_model(const Model* model, double ra, std::string filename, std::vector<double> rv);

/** This routine can be used to test and validate the implementation of new models (subclasses of the DensityModel or SurfaceDensityModel classes). It reads in a radius \f$r\f$ and calculates the density and its derivatives, the mass, and the potential at that radius. These values can be checked against values calculated in other ways. The routine also calculates the density by integrating the isotropic and the Osipkov-Merritt distribution function (with anisotropy radius \f$r_{\text{a}}\f$) over velocity space. */
void validate_model(const Model* model, double r, double ra);

/** This routine calculates the mass and the different energies for a model, that is the total potential energy, the total kinetic energy, and the total integrated binding energy. This routine also serves for validation purposes. */
void calculate_energy_model(const Model* model, double ra);

//////////////////////////////////////////////////////////////////////

#endif
