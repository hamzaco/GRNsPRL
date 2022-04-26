# Gene-Regulatory-Networks
There are two types of documents. "json" files are Python dictionaries. "py" files are Python codes. 

!!You may not able to generate all plots with given codes!!

This project includes two main class structure.

Gene Class includes:\
Attributes:\
-a boolean function\
-sensitivity value\
-activation\
-set of influencer genes\
-set of influencer externals\
-set of influencing nodes\
-set of incoming positive edges\
-set of incoming negative edges\


Functions:\
-Updating the gene:\
returns new activation\
-Replace an input with another one:\
returns updated version of genes\
-Find sensitivity:\
returns sensitivity value and number of inputs\
-Is canalizing: returns a boolean\
-Behaviour analysis betwen inputs(and-like, or-like, xor-like):\
returns number of xor-like behaviour and a boolean which detects conflicting behaviours\  
-Find symmetries on its boolean function:\
returns symmetry distrubition and a list of equivalence classes\
-Find number of truths in full logical table: return none\
-Remove an input which has no effect indeed: return none\
-Find the boolean function from its logical table and symmetries:\
returns a boolean function, and a relation which depicts the operation between classes\


Gene Network Class includes:\
Attributes:\
-expressions #txt file\
-list of externals # a list externals' names\
-list of genes\
-list of attractors with their basin size\
-coherence value\
-sensitivity value\
-correspondingExternals condition for attractors

Functions:\
-Finding Attractors: returns a list of attractors\
-Reading data from txt files\
-Initializing network relations\
-Finding equilibrium for a given initial state:\
returns the equilibrium\
-Coherence calculations\
-Sensitivity test with attractors\
-Sensitivity test with all possible states\



Logic Table Network Class -subclass of gene network- includes:\
Attributes:\
-expressions #txt file\
-list of externals # a list externals' names\
-list of genes\
-list of attractors with their basin size\
-coherence value\
-sensitivity value\
-correspondingExternals condition for attractors

Functions:\
-Finding Attractors: returns a list of attractors\
-Reading data from txt files\
-Initializing network relations\
-Finding equilibrium for a given initial state:\
returns the equilibrium\
-Coherence calculations\
-Sensitivity test with attractors\
-Sensitivity test with all possible states\



Logical Gene Class -subclass of gene network- includes:\
Attributes:\
-a logical dictionary of truth table\
-sensitivity value\
-activation\
-set of influencer genes\
-set of influencer externals\
-set of influencing nodes\
-set of incoming positive edges\
-set of incoming negative edges\

Functions:\
-Updating the gene\
-Find sensitivity\
-Construction logical dictionary from a logical function
-Is canalizing\
-Behaviour analysis betwen inputs(and-like, or-like, xor-like)
-Find symmetries on its boolean function \
-Find number of truths in full logical table\
-Find the boolean function from its logical table and symmetries\
