Installation
=================

GRASP
-------------------------------


To use GRASP, first make sure the following Matlab toolboxes are installed:

 - Parallel Computing Toolbox
 - Optimization Toolbox
 - Bioinformatics Toolbox

If you want to export models to SBML, you'll also need Matlab's SimBiology Toolbox.

If you want to use gurobi to solve the linear programs, you can get it from `here <https://www.gurobi.com/>`_.
Even though Gurobi is a commercial solver, academic licenses are free.

If you are using the ABC approach and want to get a better performance, you should also install `nlopt <https://nlopt.readthedocs.io/en/latest/>`_

Then download GRASP from `here <https://github.com/biosustain/GRASP>`_, and you're ready to go.



Python visualization
-------------------------------

If you want to use the jupyter notebooks in the visualization folder to visualize the model's Gibbs energies vs. reaction fluxes and/or the results from Metabolic Control Analysis, you'll need the following python packages:

 - jupyter notebook or jupyter lab
 - pandas
 - numpy
 - scipy
 - matplotlib
 - altair

You can also just create a conda environment and install the dependencies with poetry:

.. code-block:: shell

    conda env create -f environment.yml


then activate the conda environment with:

.. code-block:: shell

    conda activate grasp_viz


and finally install all other dependencies with poetry:

.. code-block:: shell
    
    poetry install


The commands above should be run from tho folder `visualization`.
