/*//////////////////////////////////////////////////////////////////
////     SpheCow -- Flexible dynamical models for galaxies      ////
////                and dark matter haloes                      ////
////     © Sterrenkundig Observatorium, Universiteit Gent       ////
///////////////////////////////////////////////////////////////// */

#include "SpheCow.hpp"

#include "BPLModel.hpp"
#include "DeVaucouleursModel.hpp"
#include "EinastoModel.hpp"
#include "GammaModel.hpp"
#include "GaussLegendre.hpp"
#include "HernquistModel.hpp"
#include "HypervirialModel.hpp"
#include "IsochroneModel.hpp"
#include "JaffeModel.hpp"
#include "NFWModel.hpp"
#include "NukerModel.hpp"
#include "PerfectSphereModel.hpp"
#include "PlummerModel.hpp"
#include "SersicModel.hpp"
#include "SigmoidDensityModel.hpp"
#include "SigmoidSurfaceDensityModel.hpp"
#include "ZhaoModel.hpp"

#include <iostream>
#include <fstream>
#include <iomanip>

//////////////////////////////////////////////////////////////////////

int main(void)
{
    GaussLegendre gl(128);
    double Mtot = 3.0;
    double b = 7.0;
    PlummerModel model(Mtot, b, &gl);
    double ra = 6.2;
    calculate_energy_model(&model, ra);
    return 0;
}

//////////////////////////////////////////////////////////////////////

void run_model(const Model* model, double ra, std::string filename, std::vector<double> rv)
{
    std::ofstream file(filename.c_str());
    file << "# column 0: radius" << std::endl
         << "# column 1: density" << std::endl
         << "# column 2: density slope" << std::endl
         << "# column 3: mass" << std::endl
         << "# column 4: circular velocity" << std::endl
         << "# column 5: surface density" << std::endl
         << "# column 6: surface density slope" << std::endl
         << "# column 7: surface mass" << std::endl
         << "# column 8: potential" << std::endl
         << "# column 9: isotropic dispersion" << std::endl
         << "# column 10: isotropic projected dispersion" << std::endl
         << "# column 11: isotropic distribution function" << std::endl
         << "# column 12: isotropic density of states" << std::endl
         << "# column 13: isotropic differential energy distribution" << std::endl
         << "# column 14: osipkov-merritt radial dispersion for ra = " << ra << std::endl
         << "# column 15: osipkov-merritt tangential dispersion for ra = " << ra << std::endl
         << "# column 16: osipkov-merritt projected dispersion for ra = " << ra << std::endl
         << "# column 17: osipkov-merritt distribution function for ra = " << ra << std::endl
         << "# column 18: osipkov-merritt pseudo density of states for ra = " << ra << std::endl
         << "# column 19: osipkov-merritt pseudo differential energy distribution for ra = " << ra << std::endl << std::endl;
    file << std::scientific << std::setprecision(16);
    int numr = rv.size();
    for (int i=0; i<numr; ++i)
    {
        double r = rv[i];
        std::cout << "Calculating properties for r = " << r << std::endl;
        double rho = model->density(r);
        double gamma = model->density_slope(r);
        double M = model->mass(r);
        double vc = model->circular_velocity(r);
        double Sigma = model->surface_density(r);
        double gammap = model->surface_density_slope(r);
        double Mp = model->surface_mass(r);
        double Psi = model->potential(r);
        double disp_iso = model->isotropic_dispersion(r);
        double dispp_iso = model->isotropic_projected_dispersion(r);
        double df_iso = model->isotropic_distribution_function(r);
        double g_iso = model->isotropic_density_of_states(r);
        double ded_iso = df_iso * g_iso;
        double dispr_om = model->osipkov_merritt_radial_dispersion(r,ra);
        double dispt_om = model->osipkov_merritt_tangential_dispersion(r,ra);
        double dispp_om = model->osipkov_merritt_projected_dispersion(r,ra);
        double df_om = model->osipkov_merritt_distribution_function(r,ra);
        double g_om = model->osipkov_merritt_pseudo_density_of_states(r,ra);
        double ded_om = df_om * g_om;
        file << r << '\t'
             << rho << '\t'
             << gamma << '\t'
             << M << '\t'
             << vc << '\t'
             << Sigma << '\t'
             << gammap << '\t'
             << Mp << '\t'
             << Psi << '\t'
             << disp_iso << '\t'
             << dispp_iso << '\t'
             << df_iso << '\t'
             << g_iso << '\t'
             << ded_iso << '\t'
             << dispr_om << '\t'
             << dispt_om << '\t'
             << dispp_om << '\t'
             << df_om << '\t'
             << g_om << '\t'
             << ded_om << '\t'
             << std::endl;
    }
    file.close();
    return;
}

//////////////////////////////////////////////////////////////////////

void validate_model(const Model* model, double r, double ra)
{
    std::cout << std::setprecision(12) << std::endl;
    std::cout << "rho = " << model->density(r) << std::endl;
    std::cout << "drho = " << model->derivative_density(r) << std::endl;
    std::cout << "d2rho = " << model->second_derivative_density(r) << std::endl;
    std::cout << "M = " << model->mass(r) << std::endl;
    std::cout << "Psi = " << model->potential(r) << std::endl;
    std::cout << "Mtot = " << model->total_mass() << std::endl;
    std::cout << std::endl;
    std::cout << "Density from isotropic distribution function" << std::endl;
    std::cout << "rho = " << model->density_from_isotropic_distribution_function(r) << std::endl;
    std::cout << std::endl;
    std::cout << "Density from Osikov-Merritt distribution function" << std::endl;
    std::cout << "rho = " << model->density_from_osipkov_merritt_distribution_function(r,ra) << std::endl;
    std::cout << std::endl;
    return;
}

//////////////////////////////////////////////////////////////////////

void calculate_energy_model(const Model* model, double ra)
{
    std::cout << std::setprecision(12);
    std::cout << "Properties independent of the orbital structure" << std::endl;
    std::cout << "total mass = " << model->total_mass() << std::endl;
    std::cout << "total potential energy = " << model->total_potential_energy() << std::endl << std::endl;
    std::cout << "Properties for an isotropic orbital structure" << std::endl;
    std::cout << "total mass = " << model->total_mass_from_isotropic_differential_energy_distribution() << std::endl;
    std::cout << "total kinetic energy = " << model->isotropic_total_kinetic_energy() << std::endl;
    std::cout << "total integrated binding energy = " << model->isotropic_total_integrated_binding_energy() << std::endl << std::endl;
    std::cout << "Properties for an Osipkov-Merritt orbital structure with ra = " << ra << std::endl;
    std::cout << "total mass = " << model->total_mass_from_osipkov_merritt_pseudo_differential_energy_distribution(ra) << std::endl;
    std::cout << "total kinetic energy = " << model->osipkov_merritt_total_kinetic_energy(ra) << std::endl << std::endl;
    return;
}

//////////////////////////////////////////////////////////////////////
