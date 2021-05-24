/*//////////////////////////////////////////////////////////////////
////     SpheCow -- Flexible dynamical models for galaxies      ////
////                and dark matter haloes                      ////
////     © Sterrenkundig Observatorium, Universiteit Gent       ////
///////////////////////////////////////////////////////////////// */

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

#include <iomanip>
#include <iostream>
#include <sstream>

// Python wrapper specific includes:
// Python.h (basic Python support)
#include <Python.h>
// numpy: requires the define before including the header to select an API
/*! @brief Use the NumPy 1.7 API. */
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

//////////////////////////////////////////////////////////////////////

// Helper function that converts a PyArrayObject to a 1D vector. We assume the
// array is 0D (scalar) or 1D.

std::vector<double> unpackNumpyArray(PyArrayObject *array) {
  npy_intp arraySize;
  const npy_intp arrayNdim = PyArray_NDIM(array);
  if (arrayNdim > 1) {
    throw std::string("Input array cannot be converted to scalar or 1D array!");
    // not a 0D or 1D array...
    return std::vector<double>();
  }
  if (arrayNdim > 0) {
    // actual 1D array
    const npy_intp *arrayDims = PyArray_DIMS(array);
    arraySize = arrayDims[0];
  } else {
    // the array is actually a scalar
    // we have to convert it to a 1D array
    arraySize = 1;
    npy_intp newDims[1] = {1};
    PyArray_Dims newDimsObject;
    newDimsObject.ptr = newDims;
    newDimsObject.len = 1;
    array = reinterpret_cast<PyArrayObject *>(
        PyArray_Newshape(array, &newDimsObject, NPY_ANYORDER));
  }
  std::vector<double> outputVector(arraySize);
  for (npy_intp i = 0; i < arraySize; ++i)
    outputVector[i] = *(reinterpret_cast<double *>(PyArray_GETPTR1(array, i)));

  Py_DECREF(array);
  return outputVector;
}

//////////////////////////////////////////////////////////////////////

// Helper function that converts a C++ vector to a 1D or 0D NumPy array.

PyObject *packNumpyArray(std::vector<double> &v) {
  npy_intp aDims[1] = {static_cast<npy_intp>(v.size())};
  PyArrayObject *array = reinterpret_cast<PyArrayObject *>(
      PyArray_SimpleNew(1, aDims, NPY_DOUBLE));
  for (npy_intp i = 0; i < aDims[0]; ++i)
    *reinterpret_cast<double *>(PyArray_GETPTR1(array, i)) = v[i];
  // squeezes out empty dimensions (i.e. converts a 1D to a 0D array if the
  // array only has one element)
  return PyArray_Squeeze(array);
}

//////////////////////////////////////////////////////////////////////

// get the floating point value corresponding to the given key from the given
// dictionary

double getDictElement(PyObject *dictionary, const char *key) {
  PyObject *value = PyDict_GetItemString(dictionary, key);
  if (value == nullptr) {
    std::stringstream message;
    message << "Dictionary key \"" << key << "\" not found!";
    throw message.str();
    return 0.;
  }
  double dValue = PyFloat_AsDouble(value);
  // std::cout << key << ": " << dValue << std::endl;
  return dValue;
}

//////////////////////////////////////////////////////////////////////

// most basic Pythonfunction wrapper
// note that a return value of 0 or NULL or nullptr is interpreted as an error
// by the Python interpreter not returning anything will cause the function call
// to hang in order to have no return value, you need to return None, i.e.
// Py_None.

static PyObject *run_model(PyObject *self, PyObject *args, PyObject *kwargs)
{
    const char *modelNameCString;
    PyObject *modelParameters;
    double ra;
    PyArrayObject *r;

    static char *kwlist[] = {strdup("modelName"), strdup("modelParameters"), strdup("ra"), strdup("r"), nullptr};
    // PyArray_Converter will convert any NumPy array compatible input to an array
    // object (including scalars and lists).
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "sOdO&", kwlist, &modelNameCString, &modelParameters, &ra, PyArray_Converter, &r))
    {
        // Error: wrong arguments provided
        // PyArg_ParseTupleAndKeywords will set an appropriate error code,
        // we only need to return NULL
        return nullptr;
    }

    // from here onwards we want to be able to throw exceptions and catch them
    // in order to return NULL
    try
    {
        // convert to C++ variable types
        std::string modelName(modelNameCString);
        // double parameterA = getDictElement(modelParameters, "a");
        std::vector<double> radius = unpackNumpyArray(r);

        // std::cout << "modelName: " << modelName << std::endl;
        // std::cout << "parameter: " << parameterA << std::endl;
        // std::cout << "ra: " << ra << std::endl;
        // std::cout << "radii:";
        // for (size_t i = 0; i < radius.size(); ++i) std::cout << " " << radii[i];
        // std::cout << std::endl;
        GaussLegendre *gl = new GaussLegendre(128);
        Model *model = nullptr;
        if (modelName == "BPLModel")
            model = new BPLModel(getDictElement(modelParameters, "Mtot"),
                               getDictElement(modelParameters, "rb"),
                               getDictElement(modelParameters, "beta"),
                               getDictElement(modelParameters, "gamma"), gl);
        else if (modelName == "DeVaucouleursModel")
            model = new DeVaucouleursModel(getDictElement(modelParameters, "Mtot"),
                                         getDictElement(modelParameters, "Reff"), gl);
        else if (modelName == "EinastoModel")
            model = new EinastoModel(getDictElement(modelParameters, "Mtot"),
                                   getDictElement(modelParameters, "rh"),
                                   getDictElement(modelParameters, "n"), gl);
        else if (modelName == "GammaModel")
            model = new GammaModel(getDictElement(modelParameters, "Mtot"),
                                 getDictElement(modelParameters, "b"),
                                 getDictElement(modelParameters, "gamma"), gl);
        else if (modelName == "HernquistModel")
            model = new HernquistModel(getDictElement(modelParameters, "Mtot"),
                                     getDictElement(modelParameters, "b"), gl);
        else if (modelName == "HypervirialModel")
            model = new HypervirialModel(getDictElement(modelParameters, "Mtot"),
                                       getDictElement(modelParameters, "rs"),
                                       getDictElement(modelParameters, "p"), gl);
        else if (modelName == "IsochroneModel")
            model = new IsochroneModel(getDictElement(modelParameters, "Mtot"),
                                     getDictElement(modelParameters, "b"), gl);
        else if (modelName == "JaffeModel")
            model = new JaffeModel(getDictElement(modelParameters, "Mtot"),
                                 getDictElement(modelParameters, "b"), gl);
        else if (modelName == "NFWModel")
            model = new NFWModel(getDictElement(modelParameters, "Mvir"),
                               getDictElement(modelParameters, "rs"),
                               getDictElement(modelParameters, "c"), gl);
        else if (modelName == "NukerModel")
            model = new NukerModel(getDictElement(modelParameters, "Mtot"),
                                 getDictElement(modelParameters, "Rb"),
                                 getDictElement(modelParameters, "alpha"),
                                 getDictElement(modelParameters, "beta"),
                                 getDictElement(modelParameters, "gamma"), gl);
        else if (modelName == "PerfectSphereModel")
            model = new PerfectSphereModel(getDictElement(modelParameters, "Mtot"),
                                         getDictElement(modelParameters, "c"), gl);
        else if (modelName == "PlummerModel")
            model = new PlummerModel(getDictElement(modelParameters, "Mtot"),
                                   getDictElement(modelParameters, "c"), gl);
        else if (modelName == "SersicModel")
            model = new SersicModel(getDictElement(modelParameters, "Mtot"),
                                  getDictElement(modelParameters, "Reff"),
                                  getDictElement(modelParameters, "m"), gl);
        else if (modelName == "SigmoidDensityModel")
            model = new SigmoidDensityModel(getDictElement(modelParameters, "Mtot"),
                                          getDictElement(modelParameters, "rb"),
                                          getDictElement(modelParameters, "alpha"),
                                          getDictElement(modelParameters, "beta"),
                                          getDictElement(modelParameters, "gamma"), gl);
        else if (modelName == "SigmoidSurfaceDensityModel")
            model = new SigmoidSurfaceDensityModel(getDictElement(modelParameters, "Mtot"),
                                                 getDictElement(modelParameters, "Rb"),
                                                 getDictElement(modelParameters, "alpha"),
                                                 getDictElement(modelParameters, "beta"),
                                                 getDictElement(modelParameters, "gamma"), gl);
        else if (modelName == "ZhaoModel")
            model = new ZhaoModel(getDictElement(modelParameters, "Mtot"),
                                getDictElement(modelParameters, "rb"),
                                getDictElement(modelParameters, "alpha"),
                                getDictElement(modelParameters, "beta"),
                                getDictElement(modelParameters, "gamma"), gl);
        else
        {
            std::stringstream message;
            message << "Unknown model name: \"" << modelName << "\"!";
            throw message.str();
        }

        const size_t rsize = radius.size();
        std::vector<double> density(rsize), density_slope(rsize), mass(rsize),
            circular_velocity(rsize), surface_density(rsize),
            surface_density_slope(rsize), surface_mass(rsize), potential(rsize),
            isotropic_dispersion(rsize), isotropic_projected_dispersion(rsize),
            isotropic_distribution_function(rsize),
            isotropic_density_of_states(rsize),
            isotropic_differential_energy_distribution(rsize),
            osipkov_merritt_radial_dispersion(rsize),
            osipkov_merritt_tangential_dispersion(rsize),
            osipkov_merritt_projected_dispersion(rsize),
            osipkov_merritt_distribution_function(rsize),
            osipkov_merritt_pseudo_density_of_states(rsize),
            osipkov_merritt_pseudo_differential_energy_distribution(rsize);
        for (size_t i = 0; i < rsize; ++i)
        {
            const double r = radius[i];
            density[i] = model->density(r);
            density_slope[i] = model->density_slope(r);
            mass[i] = model->mass(r);
            circular_velocity[i] = model->circular_velocity(r);
            surface_density[i] = model->surface_density(r);
            surface_density_slope[i] = model->surface_density_slope(r);
            surface_mass[i] = model->surface_mass(r);
            potential[i] = model->potential(r);
            isotropic_dispersion[i] = model->isotropic_dispersion(r);
            isotropic_projected_dispersion[i] = model->isotropic_projected_dispersion(r);
            isotropic_distribution_function[i] = model->isotropic_distribution_function(r);
            isotropic_density_of_states[i] = model->isotropic_density_of_states(r);
            isotropic_differential_energy_distribution[i] = isotropic_distribution_function[i] * isotropic_density_of_states[i];
            osipkov_merritt_radial_dispersion[i] = model->osipkov_merritt_radial_dispersion(r, ra);
            osipkov_merritt_tangential_dispersion[i] = model->osipkov_merritt_tangential_dispersion(r, ra);
            osipkov_merritt_projected_dispersion[i] = model->osipkov_merritt_projected_dispersion(r, ra);
            osipkov_merritt_distribution_function[i] = model->osipkov_merritt_distribution_function(r, ra);
            osipkov_merritt_pseudo_density_of_states[i] = model->osipkov_merritt_pseudo_density_of_states(r, ra);
            osipkov_merritt_pseudo_differential_energy_distribution[i] = osipkov_merritt_distribution_function[i] * osipkov_merritt_pseudo_density_of_states[i];
        }
        delete model;
        delete gl;
        PyObject *outputDictionary = PyDict_New();
        PyDict_SetItemString(outputDictionary, "radius", packNumpyArray(radius));
        PyDict_SetItemString(outputDictionary, "density", packNumpyArray(density));
        PyDict_SetItemString(outputDictionary, "density_slope", packNumpyArray(density_slope));
        PyDict_SetItemString(outputDictionary, "mass", packNumpyArray(mass));
        PyDict_SetItemString(outputDictionary, "circular_velocity", packNumpyArray(circular_velocity));
        PyDict_SetItemString(outputDictionary, "surface_density", packNumpyArray(surface_density));
        PyDict_SetItemString(outputDictionary, "surface_density_slope", packNumpyArray(surface_density_slope));
        PyDict_SetItemString(outputDictionary, "surface_mass", packNumpyArray(surface_mass));
        PyDict_SetItemString(outputDictionary, "potential", packNumpyArray(potential));
        PyDict_SetItemString(outputDictionary, "isotropic_dispersion", packNumpyArray(isotropic_dispersion));
        PyDict_SetItemString(outputDictionary, "isotropic_projected_dispersion", packNumpyArray(isotropic_projected_dispersion));
        PyDict_SetItemString(outputDictionary, "isotropic_distribution_function", packNumpyArray(isotropic_distribution_function));
        PyDict_SetItemString(outputDictionary, "isotropic_density_of_states", packNumpyArray(isotropic_density_of_states));
        PyDict_SetItemString(outputDictionary, "isotropic_differential_energy_distribution", packNumpyArray(isotropic_differential_energy_distribution));
        PyDict_SetItemString(outputDictionary, "osipkov_merritt_radial_dispersion", packNumpyArray(osipkov_merritt_radial_dispersion));
        PyDict_SetItemString(outputDictionary, "osipkov_merritt_tangential_dispersion", packNumpyArray(osipkov_merritt_tangential_dispersion));
        PyDict_SetItemString(outputDictionary, "osipkov_merritt_projected_dispersion", packNumpyArray(osipkov_merritt_projected_dispersion));
        PyDict_SetItemString(outputDictionary, "osipkov_merritt_distribution_function", packNumpyArray(osipkov_merritt_distribution_function));
        PyDict_SetItemString(outputDictionary, "osipkov_merritt_pseudo_density_of_states", packNumpyArray(osipkov_merritt_pseudo_density_of_states));
        PyDict_SetItemString(outputDictionary, "osipkov_merritt_pseudo_differential_energy_distribution", packNumpyArray(osipkov_merritt_pseudo_differential_energy_distribution));
        return outputDictionary;
    }
    catch (const std::string error_message)
    {
        std::cerr << error_message << std::endl;
        PyErr_SetString(PyExc_RuntimeError, error_message.c_str());
        return nullptr;
    }
}

//////////////////////////////////////////////////////////////////////

// This is a list of all the static methods that are exposed from the module.
// Each static function entry has
//  - a name (the name of the function in Python - can be the same as in C),
//  - a pointer to the corresponding C function (usually a wrapper),
//  - a type, METH_VARARGS | METH_KEYWORDS is the most general and therefore
//  recommended type,
//  - a DocString. This will show up as the __doc__ property of the function and
//  is visible from within Python.
// Note that the end of the list is (and needs to be) signaled by a row of 0s.

static PyMethodDef pySpheCowMethods[] = {
    {"run_model", reinterpret_cast<PyCFunction>(run_model),
     METH_VARARGS | METH_KEYWORDS, "Run a model."},
    {nullptr, nullptr, 0, nullptr}};

//////////////////////////////////////////////////////////////////////

// Module information struct.
// Contains all information that defines a module.
//  - The first entry needs to be set to PyModuleDef_HEAD_INIT
//  - The second entry is the name of the module. This is the name you will use
//  to import the module in Python (and the name of the library object). This
//  name needs to match the name used in the setup.py script.
//  - The third entry is the DocString for the module. This will show up as the
//  __doc__ property of the module and is visible from within Python.
//  - The fourth entry has something to do with how the module is managed in
//  memory. -1 is a sensible value here.
//  - The last entry is the list of static module functions given above.

static struct PyModuleDef pySpheCowModule = {PyModuleDef_HEAD_INIT, "pySpheCow",
                                             "SpheCow Python module.", -1,
                                             pySpheCowMethods};

//////////////////////////////////////////////////////////////////////

// Module initialisation function.
// This function is automatically called by the Python interpreter upon import
// of the module.
PyMODINIT_FUNC PyInit_pySpheCow() {
  // Initialise NumPy array support (needs to be called if you use NumPy in the
  // module).
  import_array();
  // Create the module object based on the definition given above.
  PyObject *m = PyModule_Create(&pySpheCowModule);
  // return the module object
  return m;
}

//////////////////////////////////////////////////////////////////////
