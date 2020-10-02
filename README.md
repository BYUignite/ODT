# ODT ##################

This code implements the One-Dimensional Turbulence (ODT) model for turbulent reacting or
nonreacting flows. 

## Documentation ########
Detailed documentation is available [here](https://ignite.byu.edu/ODT_documentation). 
More information on the theory and application of ODT is available
[here](https://ignite.byu.edu/ODT_documentation/odt_theory.html).

A short video overview summarizing downloading, building, running, and processing the code is shown [here]()https://vimeo.com/464356759).

<!--
The following two papers discussing theory and application of the code are available. Additional papers are available [here](http://ignite.byu.edu/publications.html).
   * [D. Lignell et al., One-dimensioanl turbulence modeling for cylindrical and spherical flows: model formulation and application, Theoretical and Computational Fluid Dynamics, 32:495-520](https://ignite.byu.edu/public/Lignell_2018.pdf)
   * [D. Lignell et al., Mesh adaption for efficient multiscale implementation of one-dimensional turbulence, Theoretical and Computational Fluid Dynamics, 27:273-295 (2013)](https://ignite.byu.edu/public/ODTmethod.pdf)
-->

## Dependencies #################

### ODT Code
* [Cantera](http://cantera.org): open-source suite of tools for problems involving chemical kinetics, thermodynamics, and transport.
* Yaml: input file format. This installation is conveniently built into the ODT build process. 
* Cmake 3.12 or higher
* (OPTIONAL) Doxygen: builds documentation. 

### Post-processing #############
Post-processing data produced by ODT and ODT is processed via Python 3 scripts. We recommend Python 3.2 or higher. Scripts may not function properly using Python 2.x. The following packages are required and can be installed via pip3:
* numpy
* scipy
* matplotlib
* glob
* yaml
* sys
* os

## Directory structure ###########
* `build`: build the code
* `data`: contains all data files output during a simulation
    * The code will generate a subfolder with a name corresponding to case name specified in the run script in the `run` folder.
        * This case subfolder will contain subfolders `input`, `runtime`, `data`, and `post`, which contain the input data files, runtime output, simulation data files, and post-processed data, respectively.
* `doc`: contains documentation files
* `input`: contains case input files
    * Other input files include a Cantera mechanism file in the `user_gas_mechanisms` folder and an optional `restart.yaml` file.
* `post`: contains post-processing scripts and files for given case types
   * Output is placed in `data/caseName/post`. These are mostly Python files. Some cases also include experimental data files for comparison and plotting.
* `run`: contains the code executable `odt.x` and several run scripts 
    * The user specifies inputDir as the path to the input file containing the case to run and specifies a case name for variable caseName. Files are created and copied into `data/caseName`, as noted above.
* `source`: contains source code (including header files) and `CMakeLists.txt` files

