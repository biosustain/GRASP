# Using GRASP with promiscuous reactions


### Defining patterns for promiscuous reactions

The flux through each reaction is generated based on the product ID.  
By product ID i mean the letters used in the pattern definition.  
For promiscuous reactions one needs to use `Pi` (`i` is an integer) as an ID for the first product to be released in each promiscuous reaction. 
`i` is related to the order of the promiscuous reactions in the stoichiometric matrix, e.g. if DDC is defined before DDC_tryptm, then the first product to be released in DDC is `P1` and in DDC_tryptm is `P2`.  



### Changes in the input excel file

The only sheet that changed was `kinetics`, where a column was added (promiscuous) to specify the sets of promiscuous reactions.  
Like inhibitors, activators, etc., the reactions are separated by a space.  
Example: for DDC the list of promiscuous is `DDC DDC_tryptm`, where the order does not matter.  
