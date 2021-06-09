/*//////////////////////////////////////////////////////////////////
////     SpheCow -- Flexible dynamical models for galaxies      ////
////                and dark matter haloes                      ////
////     © Sterrenkundig Observatorium, Universiteit Gent       ////
///////////////////////////////////////////////////////////////// */

#include "SersicModel.hpp"
#include <fstream>

//////////////////////////////////////////////////////////////////////

SersicModel::SersicModel(double Mtot, double Reff, double m, const GaussLegendre* gl)
{
    _Mtot = Mtot;
    _Reff = Reff;
    _m = m;

    // Determine the value of b. First search if from the file Sersicb.txt that
    // lists exact values for all m between 0.01 and 10 with a spacing of 0.001.
    
    std::ifstream file("Sersicb.txt");
    bool bset = false;
    double mf, bf;
    for (int j=0; j<9990; j++)
    {
        file >> mf >> bf;
        if (fabs(_m-mf)<1e-9)
        {
            _b = bf;
            bset = true;
            break;
        }
    }
    
    // If b is not yet set (it is not in the file), determine the value using the
    // approximations from Ciotti & Bertin (1999) for m>1, or from
    // Baes & Ciotti (2019) for m>1.
    
    if (!bset)
    {
        double m2 = m*m;
        if (m<1)
            _b = pow(1.0/sqrt(2.0)-0.45807*m+1.83247*m2
                         -1.2556*m2*m+0.85239*m2*m2,1.0/m);
        else
            _b = 2.0*m - 1.0/3.0 + 4.0/(405.0*m) + 46.0/(25515*m2)
                 + 131.0/(1148175*m*m2) - 2194697.0/(30690717750.0*m2*m2);
    }
    _Sigma0 = _Mtot/(_Reff*_Reff) * pow(_b,2.0*m)/(2.0*M_PI*m*tgamma(2.0*m));
    _gl = gl;
}

//////////////////////////////////////////////////////////////////////

double SersicModel::scale_radius() const
{
    return _Reff;
}

//////////////////////////////////////////////////////////////////////

double SersicModel::surface_density(double R) const
{
    double t = R/_Reff;
    double z = pow(t,1.0/_m);
    return _Sigma0 * exp(-_b*z);
}

//////////////////////////////////////////////////////////////////////

double SersicModel::derivative_surface_density(double R) const
{
    double t = R/_Reff;
    double z = pow(t,1.0/_m);
    double ef = exp(-_b*z);
    double ff = -(_Sigma0*_b) / (_m*_Reff);
    return ff * ef * (z/t);
}

//////////////////////////////////////////////////////////////////////

double SersicModel::second_derivative_surface_density(double R) const
{
    double t = R/_Reff;
    double z = pow(t,1.0/_m);
    double ef = exp(-_b*z);
    double ff = (_Sigma0*_b) / pow(_m*_Reff,2);
    return ff * ef * (z/t/t) * (-1.0+_m+_b*z);
}

//////////////////////////////////////////////////////////////////////

double SersicModel::third_derivative_surface_density(double R) const
{
    double t = R/_Reff;
    double z = pow(t,1.0/_m);
    double ef = exp(-_b*z);
    double ff = -(_Sigma0*_b) / pow(_m*_Reff,3);
    return ff * ef * (z/t/t/t) * (1.0-3.0*_m+2.0*_m*_m + 3.0*_b*(_m-1.0)*z+_b*_b*z*z);
}

//////////////////////////////////////////////////////////////////////

double SersicModel::total_mass() const
{
    return _Mtot;
}

//////////////////////////////////////////////////////////////////////

double SersicModel::central_potential() const
{
    return _Mtot/_Reff * (2.0*pow(_b,_m)*tgamma(_m+1.0)) / (M_PI*_m*tgamma(2.0*_m));
}

//////////////////////////////////////////////////////////////////////
