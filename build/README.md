## Building ODT source code with cmake

------------------------------------------------------------------
### SOFTWARE 

Required software:
    cmake (3.12 or higher)
    cantera (visit https://cantera.org/ for information and installation instructions)
    git (for installing yaml)

Optional software:
    doxygen (for building documentation)
    pdflatex (for generating PDF documentation with doxygen)

------------------------------------------------------------------
### Build instructions 

STEP 1: run cmake
    Edit the user_config file for settings and paths.
    RUN: `cmake -C user_config ../source`

STEP 2 (optional if yaml is already installed): build and install yaml
    RUN: `make yaml`

STEP 3: build the ODT code
    RUN: make -j8

(OPTIONAL) STEP 4:  build documentation
    RUN: `make doxygen`
    
------------------------------------------------------------------

### Cleanup instructions 

Basic cleanup:
    RUN: `make clean`

Thorough cleanup:
    RUN: `./clean_this_dir.sh`

------------------------------------------------------------------

### Notes

- `CMakeLists.txt` files are located in `../source` directory and its subdirectories.
- All files in this folder can be deleted except for user_config and this README.
- Generating PDF documentation requires pdflatex to be installed on your system. 

