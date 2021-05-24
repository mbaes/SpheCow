# install script for the Example module
# to build the module and place it within the current folder, run
#   python(3) setup.py build --build-lib .
# Intermediate files will be stored in a newly created build/ folder
# The module itself is placed within the current folder and can be imported in Python using
#   import Example

# imports: setuptools are required for the compilation of the module
from setuptools import setup, Extension

# glob is used to automatically list source files
import glob

# numpy is required to find numpy headers
import numpy

# main function executed when calling the script
if __name__ == "__main__":
    # list all cpp and hpp files
    cpps = glob.glob("*.cpp")
    hpps = glob.glob("*.hpp")
    # run the setup command
    # the name is not used as far as I know
    setup(
        name="PySpheCow",
        # since we need to compile our C++ extension, we add it as an extension object
        ext_modules=[
            Extension(
                "pySpheCow",  # Name of the extension. Should match the module name in the C++ code.
                cpps,  # .cpp source files that need to be compiled and linked
                include_dirs=[
                    numpy.get_include()
                ],  # additional include directories needed during compilation
                depends=hpps,  # extra dependencies for the build. By listing the header files, setuptools will recompile the module if only a header file changed.
                extra_compile_args=[
                    "-std=c++14"
                ],  # additional compilation flags passed on to the compiler.
            )
        ],
        install_requires=[
            "numpy"
        ],  # list module dependencies that need to be satisfied on the Python end
    )
