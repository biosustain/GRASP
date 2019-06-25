# GRASP

## Introduction 

GRASP is a matlab package to sample thermodynamically feasible kinetic model ensembles.
As input it takes an excel file containing the model stoichiometry, a flux distribution, and Gibbs energies.

It also performs Metabolic Control Analysis (MCA) and stability analysis on the resulting ensembles. The resulting ensembles can also be simulated. 

## Installation

### Matlab requirements

* Matlab's Parallel Computing Toolbox
* Matlab's Optimization Toobox
* Matlab's Bioinformatics Toolbox
* Gurobi

### Python requirements (only for the jupyter notebooks in `python_code`, which are not really part of GRASP)

* numpy
* pandas
* scipy
* altair (only for the jupyter notebook)
* vega (dependency of altair)
* [set_up_grasp_models](https://github.com/biosustain/set_up_grasp_models)


## Usage

### Specifying the input
To build a model ensemble using GRASP you first need to have an excel input file and respective pattern files.


#### How to specify the input and pitfalls

* The order of metabolites in every sheet of the input excel file needs to be the same as in the stoich sheet.
* The order of reactions in every sheet of the input excel file needs to be the same as in the stoich sheet.
* When setting the order of reactions in every sheet the order of reactions should be:
	1. all reactions with enzymatic mechanisms (Uni-Uni, orderedBiBi, etc.);
	2. massAction and diffusion;
	3. all reactions with no kinetic parameters in the order:
		1. freeExchange;
		2. fixedExchange.
* In the kinetics sheet
	* if there are more than one inhibitor/activator/effector these should be separated by a single space;
	* substrate and product order columns should be filled in for all reactions with enzymatic mechanisms.
* When specifying pattern files, the convention is:
	* A, B, C, D refer to substrates;
	* P, Q, R, S refer to products;
	* I refers to inhibitors, if there is more than one, try I1, I2, I3 and make sure it worked by looking into the reactions functions;
	* Z refers to the activators, if there is more than one, try Z1, Z2, Z3 and make sure it worked by looking into the reactions functions.
* If setting compute thermodynamics to 1, then in thermoRxns specify the Gibbs standard energies, otherwise specify the reaction Gibbs energies.
* When including isoenzymes, in the rxns tab call each of the isoenzymes by their 'parent' enzyme name e.g. 'HEX' for HEX1 and HEX2:
	* they will be grouped by this name so make sure none are conflicting;
	* it is not essential to include names if there aren't isoenzymes, however the isoenzyme column is neccessary;
	* this column can also be used for promiscuous reactions if you don't know how much of the total flux goes through each reaction.


[set_up_grasp_models](https://github.com/biosustain/set_up_grasp_models) can also be used to generate pattern files.


#### Checking your GRASP input excel file

Before running GRASP it is a good idea to use [set_up_grasp_models](https://github.com/biosustain/set_up_grasp_models) to check your GRASP input file.



#### How to make sure reaction mechanisms are specified correctly:

1. go into the folder reactions_<modelID>, find the `model_name.m` file and open it;
2. check the order of metabolites in the flux function call;
3. open the respective flux function and make sure that A,B,P,Q, etc. correspond to the metabolites you expect.


#### Defining patterns for promiscuous reactions

The flux through each reaction is generated based on the product ID. By product ID i mean the letters used in the pattern definition.  

For promiscuous reactions one needs to use `Pi` (`i` is an integer) as an ID for the first product to be released in each promiscuous reaction.
`i` is related to the order of the promiscuous reactions in the stoichiometric matrix, e.g. if DDC is defined before DDC_tryptm, then the first product to be released in DDC is `P1` and in DDC_tryptm is `P2`.  


## Building models

Once the input excel file and respective pattern files are done, run `build_model.m`.

Don't forget to change the modelID to be the name of your input excel file, as well as the input and output folders in `inputFile` and `outputFile`.

In principle now GRASP only outputs models that are stable and thermodynamically feasible. This means that models which don't satisfy all the constraints are discarded and GRASP will keep sampling models until the number of (stable) models specified in the excel file is reached. 
To avoid that GRASP keeps sampling models forever (can happen when it doesn't sample any stable model), a new variable `maxNumberOfSamples` was introduced to set an upper limit to how many models can actually be sampled, irrespective of how many are stable. 

## Analyzing model ensembles

Once the model ensemble is generated you might want to analyze the results:
 - to do Metabolic Control Analysis (MCA), use the `MCA_analysis.m` script;
 - to simulate your model use the `simulate_model.m` script;
 - to check the stability of your model, i.e. the eigenvalues of the jacobian, use the `stability_analysis.m`;


**MCA analysis**

To do MCA analysis on your model ensemble, all you need to do is change the `modelID` and `outputFolder` variables in `MCA_analysis.m` and run the script.

If you set `saveMCAMatrices` to 1, all the MCA results for each model in the ensemble will be stored in the `mcaResults` variable, otherwise you will only get average values for the whole ensemble.
Note that for large models saving the results for each model will increase the runtime and require a significant amount of free space in your hard drive.

If you have promiscuous enzymes you might want to use the function `controlAndResponseAnalysis` instead of `controlAnalysis`, since response analysis answers the question: "when increasing the concentration of a given enzyme how does the steady-state flux changes?", while control analysis answers the question "when increasing the flux through a given reaction how does the steady-state flux changes?"

To visualize the results you might want to use the jupyter notebook `visualize_controlAnalysis_matlab` in the `python_code` folder.


**Model simulation**

To simulate your model you first need to convert the matlab file in which the model is encoded into a format that can be used by the Matlab ODE solvers. 
That model file is usually found in the reactions folder and is named as `modelID_Kinetics1.m`

To do so you can either do it manually or use [set_up_grasp_models](https://github.com/biosustain/set_up_grasp_models)'s function `convert_to_ode_model` as follows: 

```python
from set_up_grasp_models.set_up_models.set_up_ode_model import convert_to_ode_model

file_in = 'reactions_model_v2_3_all/model_v2_3_all_Kinetics1.m'
convert_to_ode_model(file_in)
```

It will create a new file: `reactions_model_v2_3_all/model_v2_3_all_Kinetics1_ode.m`.

The `simulateEnsemble` functions expects to find the ODE file in the reactions folder of the respective model (named as `reactions_modelID`) with the name `modelID_Kinetics1_ode.m`.


## Know issues and limitations

 - massAction mechanisms only work for one substrate and one product;
 - if the Gibbs energy range for a given reaction is large enough, it is possible that GRASP samples Gibbs energies incompatible with the reference flux. These models are discarded at the moment;
 - when using altair in the jupyter notebooks you might get the error `<VegaLite 2 object>` when trying to plot something. 
    - you should add `alt.renderers.enable('default')` after importing altair. 
    
    
 


