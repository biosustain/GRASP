# GRASP 

Branch with support for promiscuous reaction.  
Note that this branch is still experimental and is not guaranteed to work for the example included.  
Use at your own risk.  


## Using GRASP with promiscuous reactions


### Defining patterns for promiscuous reactions

The flux through each reaction is generated based on the product ID.  
By product ID i mean the letters used in the pattern definition.  
For promiscuous reactions one needs to use `Pi` (`i` is an integer) as an ID for the first product to be released in each promiscuous reaction. 
`i` is related to the order of the promiscuous reactions in the stoichiometric matrix, e.g. if DDC is defined before DDC_tryptm, then the first product to be released in DDC is `P1` and in DDC_tryptm is `P2`.  



### Changes in the input excel file

The only sheet that changed was `kinetics`, where a column was added (promiscuous) to specify the sets of promiscuous reactions.  
Like inhibitors, activators, etc., the reactions are separated by a space.  
Example: for DDC the list of promiscuous is `DDC DDC_tryptm`, where the order does not matter.  




## General overview of code changes


**Gibbs energies**

* are now stored in `ensemble` so that they can be re-used to calculate the reversibilities


**Reversibilities**

* if reactions share any intermediate steps, it's the same as for random mechanisms;
* if reactions do not share any intermediate steps, the reversibilities must sum to 1 for each reaction;
* for both, make sure that each set of promiscuous reactions has the same reversibilities.


**Allostery**

* reactionFlux now stored in ensemble.reactionFluxAllosteric
* so far sampling is done independently for each reaction, promiscuous or not


*Gibbs energies, reversibilities, and allostery are now calculated before the main loop. This is because when doing certain calculations for promiscuous reaction i one needs to update/get values from reaction j, where j > i.*


**Enzyme intermediate abundances**

* if it's a promiscuous reaction and it' not the first one in the set of promiscuous reactions, just get the abundances from the first one


**Branch factor**

* if the reaction is promiscuous set branchFactor to the promiscuous reactions fluxes and reactionFlux as the sum of all promiscuous reactions fluxes
* when having promiscuous reaction with inhibitions, the inhibition entry in reverTemp is 0 and not 1 as for other reactions (accounting for that)


**fluxReaction**

* if the reactions share common steps, the fluxReaction should be the sum of the promiscuous reactions fluxes
* if the reactions are independent, the fluxReaction should be the maximum flux of the promisucous reactions
* this is because otherwise the predicted fluxes will be wrong.


**cycledPaths**

* now it also stores paths depending on the prod node (`paths{i, 1} -> paths{i, j}`)
* to avoid weird paths on AANAT when assuming independent promiscuous reactions, added condition `subsNodes(i) < prodNodes(j)`


*Promiscuous reactions that share intermediate steps are actually just a special case of random mechanisms, assuming the gibbs energy is the same*


**Notes from 27/6/2018**

Problem: kinetic parameters for the promiscuous reactions are different (e.g. for DDC/AANAT and DDC/AANAT_tryptm) given the same elemenFlux and enzymeVect. This is because the flux is not the same for these promiscuous reactions.

Option 1: when calculating the kinetic parameters set the flux to the sum of both fluxes, e.g. v_DDC + v_DDC_tryptm
 - issue: one loses information about the fluxes going through each branch.

Option 2:  calculate the kinetic parameters for each branch individually using the respective flux.
  - issue: when there are common rate constants for the two reactions it probably won't work.

Option 3: define the branching factor based on the flux through the reactions, instead of sampling it and define the flux through the reaction as the sum of both fluxes, e.g. v_DDC + v_DDC+tryptm
