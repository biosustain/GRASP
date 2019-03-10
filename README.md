# GRASP

## Matlab requirements

* Matlab's Parallel Computing Toolbox
* Matlab's Optimization Toobox
* Matlab's Bioinformatics Toolbox
* Gurobi

## Python requirements (only for the scripts in `python_code`, which are not really part of GRASP)

* numpy
* pandas
* altair (only for the jupyter notebook)
* vega (dependency of altair)

## Checking your GRASP input Excel file

To make sure your GRASP input file has no issues (or at least a few less) before running GRASP, you can use:
 1. the script `check_input_model.py` in `python_code` to check:
    1. that the metaboltite/reaction order in all sheets is the same as in stoic;
    2. that, in the kinetics sheet, there are no metabolites separated by commas, semi colons, or dots (these should be separated by a space).
 2. `explore_initial_dGs_and_fluxes` in `python_code`, to calculate and visualize each reaction flux and delta G (given a GRASP input model). It is useful to find the TMFA problematic reactions.


## How to specify the input and pitfalls

* The order of metabolites in every sheet of the input excel file needs to be the same as in the stoich sheet;
* the order of reactions in every sheet of the input excel file needs to be the same as in the stoich sheet;
* when setting the order of reactions in every sheet the order of reactions should be:
	1. all reactions with kinetic parameters (Uni-Uni, massAction, etc.);
	2. all reactions with no kinetic parameters in the order:
		1. freeExchange
		2. fixedExchange
* in the kinetics1 sheet
	* if there are more than one inhibitor/activator/effector these should be separated by a single space;
	* if you don't specify the order of ligand binding/release it will be taken from the stoichiometric matrix;
* when specifying pattern files, the convention is:
	* A,B,C refer to substrates
	* P,Q,R refer to products
	* I refers to inhibitors, if there is more than one, try I1, I2, I3 and make sure it worked by looking into the reactions functions.
* if setting compute thermodynamics is set to 1, then in thermoRxns specify the Gibbs standard energies, otherwise specify the reaction Gibbs energies.
* when including isoenzymes, in the rxns tab call each of the isoenzymes by their 'parent' enzyme name e.g. 'HEX' for HEX1 and HEX2.
	* They will be grouped by this name so make sure none are conflicting


To check if the order of metabolites and reactions in all sheets is consistent with stoic run your model through `check_input_model.py`. This script will also check if there are any commas, semi-colons, and dots on the metabolites list in the kinetics sheet.


### How to make sure reaction mechanisms are specified correctly:

1. Go into the folder reactions1, find the `model_name.m` file and open it;
2. check the order of metabolites in the flux function call;
3. open the respective flux function and make sure that A,B,P,Q, etc correspond to the metabolites you expect.


## Using GRASP with promiscuous reactions


### Defining patterns for promiscuous reactions

The flux through each reaction is generated based on the product ID.  
By product ID i mean the letters used in the pattern definition.  
For promiscuous reactions one needs to use `Pi` (`i` is an integer) as an ID for the first product to be released in each promiscuous reaction.
`i` is related to the order of the promiscuous reactions in the stoichiometric matrix, e.g. if DDC is defined before DDC_tryptm, then the first product to be released in DDC is `P1` and in DDC_tryptm is `P2`.  
