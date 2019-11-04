.. GRASP documentation master file, created by
   sphinx-quickstart on Mon Sep 30 18:05:30 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to GRASP's documentation!
=================================


GRASP is a framework used to build kinetic model ensembles based on a single data point, which is characterized by a flux distribution and reaction free Gibbs energies.

The resulting model ensemble can then be simulated or subject to Metabolic Control Analysis.

What do you actually need to build the model?

 - fluxes for all reactions;
 - standard Gibbs energies for all reactions;
 - absolute metabolite concentrations;
 - if available, enzyme regulation information (enzyme inhibitors and activators);
 - if available, enzyme mechanism information.


For more information on how the framework works, please see

 - https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004195
 - https://www.nature.com/articles/srep29635



Documentation contents
-----------------------

.. toctree::
   :maxdepth: 1

   installation
   how_to_build_a_model
   model_analysis
   visualization
   tutorials
   api_doc
   how_to_contribute

