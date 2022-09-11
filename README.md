# GRASP

<a href='https://graspk.readthedocs.io/en/latest/?badge=latest'>
    <img src='https://readthedocs.org/projects/graspk/badge/?version=latest' alt='Documentation Status' />
</a>

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
- [Create visualizations environment](#create-visualizations-environment)
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


## Create visualizations environment

To use the jupyter notebooks for result visualization (under `visualization` folder), make sure you first go to the folder `visualization`, create the conda environment and install the necessary dependencies with:

```
conda env create -f environment.yml
```

then activate the conda environment with:

```
conda activate grasp_viz
```

and finally install all other dependencies with poetry:

```
poetry install
```


## Building the documentation

The documentation is done in sphinx, to build it locally start by going to the folder `docs` and creating a conda environment with:

```
conda env create -f environment.yml
```

then activate the conda environment with:

```
conda activate grasp_docs
```

and finally install all other dependencies with poetry:

```
poetry install
```

Once all the dependencies are installed, you can just run `make html`  while on the `docs` folder and that will create a folder `_build/html` where you can find the html documentation.


## Running unit tests

To run the unit tests go to `matlab_code/tests` and run the script `run_unit_tests.m`.


## Know issues and limitations

 - massAction mechanisms only work for one substrate and one product;
 - if the uncertainty in the flux is such that it can be both positive and negative, it is possible that GRASP samples Gibbs energies incompatible with the reference flux. This generates an error;
 - when using altair in the jupyter notebooks you might get the error `<VegaLite 2 object>` when trying to plot something. 
    - you should add `alt.renderers.enable('default')` after importing altair. 
