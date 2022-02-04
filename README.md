# ODT ##################

This code implements the One-Dimensional Turbulence (ODT) model for turbulent reacting or nonreacting flows. See also the [Basic ODT](https://github.com/BYUignite/basicODT) implementation.

## Documentation ########
Detailed documentation is available [here](https://ignite.byu.edu/ODT_documentation). 
More information on the theory and application of ODT is available
[here](https://ignite.byu.edu/ODT_documentation/odt_theory.html).

A short video overview summarizing downloading, building, running, and processing the code is shown [here](https://youtu.be/unsMJiDpSVY).

<!--
The following two papers discussing theory and application of the code are available. Additional papers are available [here](http://ignite.byu.edu/publications.html).
   * [D. Lignell et al., One-dimensioanl turbulence modeling for cylindrical and spherical flows: model formulation and application, Theoretical and Computational Fluid Dynamics, 32:495-520](https://ignite.byu.edu/public/Lignell_2018.pdf)
   * [D. Lignell et al., Mesh adaption for efficient multiscale implementation of one-dimensional turbulence, Theoretical and Computational Fluid Dynamics, 27:273-295 (2013)](https://ignite.byu.edu/public/ODTmethod.pdf)
-->

## Dependencies #################

### ODT Code
* Cmake 3.14+
* SCons 3.3.0+
* Boost 1.55+
* Git 2.0+
* (OPTIONAL) Doxygen 1.8+

### Post-processing Tools #############
Post-processing data produced by ODT and ODT is processed via Python 3 scripts. We recommend Python 3.2 or higher. Scripts may not function properly using Python 2.x. The following packages are required and can be installed via pip3:
* numpy
* scipy
* matplotlib
* glob
* yaml
* sys
* os


## Building ODT

Prior to building the ODT code, ensure that you have the required [dependencies](@ref dependencies) installed on your machine. These instructions assume that you are building ODT from a Linux-like command line.

1. Create and navigate into a `build` directory.
2. Configure CMake: `cmake ..`
3. Build and install Cantera: `make cantera`
4. Build the ODT code: `make`
5. Install ODT: `make install`

### Notes

- CMake options can be specified in the `CMakeCache.txt` file or on the command line as follows: `cmake .. -D[ODT_OPTION]=ON/OFF`
- Cantera does not need to be rebuilt if changes are made to the CMake configuration or the ODT code. 