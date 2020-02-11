# SEC (Stochastic Eddy Cascade)

This code implements the One-Dimensional Turbulence (ODT) model for turbulent reacting or nonreacting flows. 

## Documentation
Detailed documentation is available [here](http://ignite.byu.edu/SEC). The following two papers discussing theory and application of the code are available. Additional papers are available [here](http://ignite.byu.edu/publications.html).
   * [D. Lignell et al., One-dimensioanl turbulence modeling for cylindrical and spherical flows: model formulation and application, Theoretical and Computational Fluid Dynamics, 32:495-520](https://ignite.byu.edu/public/Lignell_2018.pdf)
   * [D. Lignell et al., Mesh adaption for efficient multiscale implementation of one-dimensional turbulence, Theoretical and Computational Fluid Dynamics, 27:273-295 (2013)](https://ignite.byu.edu/public/ODTmethod.pdf)

## Directory structure
* `build_make`: build the code using standard GNU make.
* `cmake_build`: build the code using cmake.
* `data`: contains all data files output during a simulation.
    * The code will generate a subfolder with a name corresponding to case name specified in the run script in the run folder.
        * This case subfolder will contain subfolders `input`, `runtime`, `data`, and `post`, which contain the input data files, runtime output, simulation data files, and post-processed data, respectively.
* `doc`: contains documentation files built using Doxygen.
* `input`: contains case input files, notably `input.yaml`, which is the primary input file with code parameters.
    * Other input files are a Cantera mechanism file in `gas_mechanisms` and an optional `restart.yaml` file.
* `post`: contains post-processing scripts and files for given cases. 
   * Output is placed in `data/caseName/post`. These are mostly Python files. Some cases also include experimental data files for comparison and plotting.
* `run`: contains the code executable `sec.x` and several run scripts. These scripts are run by the user to execute the code.
    * The user specifies inputDir as the path to the input file containing the case to run and specifies a case name for variable caseName. Files are created and copied into `data/caseName`, as noted above.
    * The user chooses one of the following run scripts to execute the code: 
      * `runOneRlz.sh` will run a single realization of the code. This is appropriate for some cases, like a statistically stationary channel flow. Many cases require running many realizations to gather turbulent statistics.
      * `runManRlz.sh` will run many realizations on a single processor.
      * `slrmJob.sh` is an example script for running embarrasingly parallel simulations, i.e. one realization for each MPI process.
      * `slrmJob_array.sh` is an example script that runs multiple realizations on a parallel machine using a slurm array.
    * `changeInputParam.py` is a convenience script for changing a value of a variable in the input file. This is convenient when running several cases changing an input parameter and can be used within the run scripts listed above.
* `source`: contains source code (including header files) and `CMakeLists.txt` files.

## Dependencies
### SEC code
* [Cantera](http://cantera.org): open-source suite of tools for problems involving chemical kinetics, thermodynamics, and transport.
* Yaml: input file format. This installation is conveniently built into the SEC build process. 
* Cmake 3.12 or higher
* (OPTIONAL) Doxygen: builds documentation. 
### Post-processing
Post-processing data produced by SEC and ODT is processed via Python 3 scripts. We recommend Python 3.2 or higher. Scripts may not function properly using Python 2.x. The following packages are required and can be installed via pip3:
* numpy
* scipy
* matplotlib
* glob
* yaml
* sys
* os

## Build
Two build systems are available: a standard GNU make, and cmake. We recommend the cmake version whenever possible. See the README files in the `build_make` and `cmake_build` folders for details.

## Test cases
### Channel flow
  1. Build the code using either Cmake or GNU make. We recommend Cmake whenever possible.
  2. Navigate to the `run` directory. Open `runOneRlz.sh` and confirm that the input file path and case name are set properly. The defaults in `runOneRlz.sh` are `inputDir="../input/channelFlow"` and `caseName="channel"`.
  3. Run `./runOneRlz.sh` to run the case. It should take less than two minutes on an average system. 
  4. Navigate to `post/channelFlow`. 
  5. Run `python3 stats.py [caseName]`. With the default case name, this becomes `python3 stats.py channel`. This will generate two plots in `../data/[caseName]/post`. Navigate there to view them. 
  6. In `../data/[caseName]/post`, there should be two newly-generated PDFs. These two plots compare the mean and RMS velocity profiles of the channel flow case just run with SEC to previous DNS data of the same case. 
### Reacting jet
  1. As before, build the code using either Cmake or GNU make. In `user_config`, change the `CHEMISTRY` parameter to `SIMPLEDLR`. If the code was previously built with this parameter, it does not need to be rebuilt. 
  2. Navigate to the `run` directory. 
      * To run one realization, open `runOneRlz.sh` and change the parameters to `inputDir="../input/jetFlame/DLR_A"` and `caseName="reacting_jet"`. Then run `./runOneRlz.sh` to start the case. On average, this should take 15-20 minutes. The resulting data can be visualized, but it is not very useful on its own. Since ODT is a stochastic model, it works best when many realizations are averaged together. 
      * To run multiple realizations, instead open `runManyRlzs.sh` and change the parameters to `nRlz=8`, `inputDir="../input/jetFlame/DLR_A"`, and `caseName="reacting_jet"`. Run `./runManyRlz.sh` to run the case. This will run eight individual realizations in series, where each will take 15-20 minutes to complete. Since we typically run several hundred realizations per simulation, this is more efficiently done in parallel using supercomputing resources, but small numbers of realizations are still instructive for test cases. 
  3. Navigate to `post/jetFlame/DLR_A`. 
  4. Run `python3 driver.py [caseName]`. For this case, the command becomes `python3 driver.py reacting_jet`. This will generate various data files and plots in `../data/[caseName]/post`. Navigate there to view them. 
  5. In `../data/[caseName]/post`, there should be several newly-generated PDFs and text data files. For example, `cl_uvel_reacting_jet.pdf` will show the axial velocity and its RMS along the centerline at various heights in the jet, `cl_temp_reacting_jet.pdf`is the corresponding centerline gas temperature, and `cl_mixf_reacting_jet.pdf`the centerline mixture fraction. 
