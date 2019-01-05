# Equation-free modeling
In this repository, I work on developing a module that contains code from the equation-free modeling developed by Kevrekidis et al. (2004) and refined by Marschler et al. (2014). Up until now, my focus was on developing a method for a bifurcation analysis using the equation-free method (as described in Marschler et al. (2014)). In the example directory, examples for the application of the equation-free modeling are given.

NOTE: This is work in progress. But the current version is working for the examples given in the examples folder. I am very happy about any feedback and contribution to the project!

# General concept
Writing the code for the equation-free methods (e.g. bifurcation analysis in Marschler et al. (2014)) can be a tedious job. In this project, I want to develop a framework such that equation-free methods are easily applicable to any microscopic model (e.g. traffic model, public goods games, etc). The eventual aim is, that a user just needs to provide the lifting, evolution and restriction operator and the parameter settings of the considered microscopic model to make use of methods such as bifurcation analysis or projective integration. 

Currently, this is done as follows: The user defines for each operator (lifting, evolution, restriction) functions that perform the particular task. The lifting operator turns a macroscopic state into a microscopic state (function takes macroscopic state dictionary and returns microscopic state as python dictionary). The evolution operator changes the microscopic state (function takes microscopic state dictionary and returns microscopic state as python dictionary). The restriction operator turns a microscopic state into a macroscopic state (function takes microscopic state dictionary and returns macroscopic state as python dictionary). The function is written as if they were already part of the eqfModel instance. Hence, the first argument is "self". This opens up the possibility to make use of a previously computed reference state (can be used i.e. in the lifting operator). Python dictionaries are used because they make it easy to work in the in the eqfModel object in a dynamic way with the state variables.

# Literature
Kevrekidis, I. G., Gear, C. W., & Hummer, G. (2004). Equation‐free: The computer‐aided analysis of complex multiscale systems. AIChE Journal, 50(7), 1346-1355.

Marschler, C., Sieber, J., Berkemer, R., Kawamoto, A., & Starke, J. (2014). Implicit methods for equation-free analysis: convergence results and analysis of emergent waves in microscopic traffic models. SIAM Journal on Applied Dynamical Systems, 13(3), 1202-1238.
