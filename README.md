# GRASP

## Introduction 

GRASP is a Matlab package to sample thermodynamically feasible kinetic model ensembles.
As input it takes an excel file containing the model stoichiometry, a flux distribution, and Gibbs energies.

It also performs Metabolic Control Analysis (MCA) and stability analysis on the resulting ensembles. The resulting ensembles can also be simulated, and individual models can be exported to SBML.

The documentation can be found [here](https://graspk.readthedocs.io).


## Table of contents

- [Introduction](#introduction)
- [Installation](#installation)
  - [Matlab requirements](#matlab-requirements)
  - [Python requirements (only for the jupyter notebooks in `visualization`, which are not really part of GRASP)](#python-requirements-only-for-the-jupyter-notebooks-in-visualization-which-are-not-really-part-of-grasp)
- [Usage](#usage)
- [Building the documentation](#building-the-documentation)
- [Running unit tests](#running-unit-tests)
- [Know issues and limitations](#know-issues-and-limitations)


## Installation

### Matlab requirements

* Matlab's Parallel Computing Toolbox
* Matlab's Bioinformatics Toolbox
* Matlab's Optimization Toobox (to run GRASP in parallel mode)
* Matlab's SimBiology Toolbox (to export models to SBML)

If you don't want to use the Matlab solver to solve the linear programming problems, you can install Gurobi and just specify in the input file (`general` sheet) that you want to use gurobi instead of linprog.

### Python requirements (only for the jupyter notebooks in `visualization`, which are not really part of GRASP)

* jupyter lab or jupyter notebook
* numpy
* pandas
* scipy
* matplotlib
* altair

## Usage

GRASP allows you to:
 - build a model ensemble;
 - do metabolic control analysis on the ensemble;
 - simulate the ensemble.

To do these, you can find example scripts in the `examples` folder and tutorials in the `tutorials` folder.
 
Once you have the results, you can use the jupyter notebooks in the `visualization` folder to visualize the results.

For more details, see the documentation [here](https://graspk.readthedocs.io).


## Building the documentation

The documentation is done in sphinx, so first you need to install all the dependencies listed in the requirements file `requirements.txt` in `docs`. 
Once all the dependencies are installed (please make sure the versions are correct, very important!), you can just run `make html`  while on the `docs` folder.


## Running unit tests

To run the unit tests go to `matlab_code/tests` and run the script `run_unit_tests.m`.


## Know issues and limitations

 - massAction mechanisms only work for one substrate and one product;
 - if the uncertainty in the flux is such that it can be both positive and negative, it is possible that GRASP samples Gibbs energies incompatible with the reference flux. This generates an error;
 - when using altair in the jupyter notebooks you might get the error `<VegaLite 2 object>` when trying to plot something. 
    - you should add `alt.renderers.enable('default')` after importing altair. 
