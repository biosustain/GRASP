.. _visualization_python:

Visualizing results
=============================================

Three jupyter notebooks are provided in the ``visualization`` folder to visualize:

 - the reactions Gibbs free energies and fluxes, to check if these are compatible;
 - the metabolic control analysis (MCA) results;
 - the results from model simulations.

The package requirements are specified in the ``requirements.txt`` file. |br|
In particular you'll need to install the following packages:

 - jupyter notebook or jupyter lab 
 - pandas
 - numpy
 - scipy
 - matplotlib
 - altair

To install the visualization package itself, open the terminal, go to the ``visualization/simulation_viz`` folder and type ``pip install .``



Visualizing fluxes and Gibbs energies
---------------------------------------

The jupyter notebook ``visualize_initial_dgs_and_fluxes`` uses `altair <https://altair-viz.github.io/>`_ to plot the Gibbs free energy and respective flux for each reaction side by side.
The goal is to easily visualize if all Gibbs energies and fluxes in the model are compatible.

The main input is the input excel file that GRASP also takes as input to build the model ensemble.

To run the notebook, all you need is to define the path to the input excel file in the section ``Define path to input excel file`` and then run all the cells.



Visualizing MCA results
---------------------------------------

The jupyter notebook ``visualize_mca`` uses `altair <https://altair-viz.github.io/>`_ to plot the control coefficients obtained from metabolic control analysis.

It uses heatmaps to plot both the median flux and concentration control coefficients along with the respective interquartile ranges.

To run the notebook, you must first specify:

 - the path to the folder that contains the MCA results;
 - the model ID (the same you defined when you built the model and ran the MCA in Matlab);
 - the number of models in the ensemble.

Once that is done, just run all the cells.



Visualizing model simulations
---------------------------------------

The jupyter notebook ``visualize_simulations`` uses `altair <https://altair-viz.github.io/>`_ and `matplotlib <https://matplotlib.org/>`_ to visualize the results of simulating GRASP models. You can choose whether to use only altair or matplotlib.

To run the notebook you need to specify:

 - the path to the folder that contains both the model ensemble ``.mat`` file and the simulation results ``.mat`` file;
 - the model name (file name for the model ensemble ``.mat`` file);
 - the simulation name (file name for simulation ``.mat`` file, without the ``simulation_`` and ``.mat`` part);
 - the number of models to be simulated.

For more info about each function in the notebook, what it does and what are the arguments, please see :ref:`simulation_viz_python`.


.. |br| raw:: html

      <br>