# GRASP 

## Requirements

* Matlab's Parallel Computing Toolbox
* Matlab's Optimization Toobox
* Matlab's Bioinformatics Toolbox
* Gurobi


## How to specify the input and pitfalls

* The order of metabolites in every sheet of the input excel file needs to be the same as in the stoich sheet;
* the order of reactions in every sheet of the input excel file needs to be the same as in the stoich sheet;
* in the kinetics1 sheet
	* if there are more than one inhibitor/activator/effector these should be separated by a single space;
	* if you don't specify the order of ligand binding/release it will be taken from the stoichiometric matrix;
* when specifying pattern files, the convention is:
	* A,B,C refer to substrates
	* P,Q,R refer to products
	* I refers to inhibitors, if there is more than one, try I1, I2, I3 and make sure it worked by looking into the reactions functions.
* if setting compute thermodynamics is set to 1, then in thermoRxns specify the Gibbs standard energies, otherwise specify the reaction Gibbs energies.


To check if the order of metabolites and reactions in all sheets is consistent with stoic run your model through `check_input_model.xlsx`. This script will also check if there are any commas, semi-colons, and dots on the metabolites list in the kinetics sheet.


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




