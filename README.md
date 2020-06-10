# GRASP

## Introduction 

GRASP is a Matlab package to sample thermodynamically feasible kinetic model ensembles.
As input it takes an excel file containing the model stoichiometry, a flux distribution, and Gibbs energies.

It also performs Metabolic Control Analysis (MCA) and stability analysis on the resulting ensembles. The resulting ensembles can also be simulated. 

For now the documentation can be accessed by opening the file `index.html` in `GRASP/docs/_build/html`  (temporary solution).


## Table of contents

* [Introduction](#introduction)
* [Installation](#installation)
* [Usage](#usage)
  * [Specifying the input](#specifying-the-input)
  * [Building models](#building-models)
  * [Analyzing model ensembles](#analyzing-model-ensembles)
* [Building the documentation](#building-the-documentation)
* [Running unit tests](#running-unit-tests)
* [Known issues and limitations](#known-issues-and-limitations)


## Installation

### Matlab requirements

* Matlab's Parallel Computing Toolbox
* Matlab's Optimization Toobox
* Matlab's Bioinformatics Toolbox
* Matlab's SimBiology Toolbox (to export models to SBML)
* Gurobi

### Python requirements (only for the jupyter notebooks in `visualization`, which are not really part of GRASP)

* numpy
* pandas
* scipy
* altair (only for the jupyter notebook)
* vega (dependency of altair)
* [set_up_grasp_models](https://github.com/biosustain/set_up_grasp_models)


## Usage

GRASP allows you to:
 - build a model ensemble;
 - do metabolic control analysis on the ensemble;
 - simulate the ensemble.

To do these, you can find example scripts in the `examples` folder.
 
Once you have the results, you can use the jupyter notebooks in the `visualization` folder to visualize the results.

For more details, see the documentation by opening the `index.html` file at `GRASP/docs/_build/html` (temporary solution).


## Building the documentation

The documentation is done in sphinx, so first you need to install all the dependencies listed in the requirements file `requirements.txt` in `docs`. 
Once all the dependencies are installed (please make sure the versions are correct, very important!), you can just run `make html`  while on the `docs` folder.


## Running unit tests

To run the unit tests, at the moment, go to `matlab_code/tests` and run the script `run_unit_tests.m`.


## Know issues and limitations

 - massAction mechanisms only work for one substrate and one product;
 - if the uncertainty in the flux is such that it can be both positive and negative, it is possible that GRASP samples Gibbs energies incompatible with the reference flux. These models are discarded at the moment;
 - when using altair in the jupyter notebooks you might get the error `<VegaLite 2 object>` when trying to plot something. 
    - you should add `alt.renderers.enable('default')` after importing altair. 
 - GRASP might no longer work with more than one model structure.
