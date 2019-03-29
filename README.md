# Equation-free modeling
In this repository, I am working on a module that contains code from the equation-free modeling developed by Kevrekidis et al. (2004). Later, this method was refined by Marschler et al. (2014). Up until now, my focus was on wirting a the code for a coarse bifurcation analysis. The theory for this is described in Marschler et al. (2014). In the example directory, I show some codes for the coarse bifuraction analysis of traffic models.

NOTE: This is work in progress. But the current version is working for the examples given in the examples folder. I am very happy about any feedback and contribution to the project!

# Requirements
NumPy

Pandas

matplotlib

# General concept
Writing the code for the equation-free modeling (e.g. bifurcation analysis in Marschler et al. (2014)) for a certain microscopic model can be a tedious job. In this project, I want to develop a framework such that equation-free modeling can easily be applied to any microscopic model (e.g. traffic model, public goods games, etc). The eventual aim is, that a user just needs to provide the lifting, evolution and restriction operator and the parameter settings of the considered microscopic model to make use of methods such as bifurcation analysis or projective integration. 

Currently, this is done as follows: The user defines for each operator (lifting, evolution, restriction) functions that perform the particular task. The function are written as if they were already part of the eqfModel instance. Hence, the first argument is "self". This opens up the possibility to make use of a previously computed reference state that is a property of the eqfModel instance. This is particularly important i.e. for the lifting operator. In general, the operator functions need such that the microscopic and macroscopic states are stored (and returned) in python dictionaries. Dictionaries are used because they make it easy to work in the in the eqfModel object in a dynamic way with the state variables.

## Lifting operator 
The lifting operator turns a macroscopic state into a microscopic state. The fucntion macroscopic state (python dictionary) and returns microscopic (python dictionary). 

## Evolution operator
The evolution operator takes the microscopic state and evolves it in time using the microscopic model. The function takes microscopic state (python dictionary) and returns microscopic state (python dictionary). 

## Restriction operator
The restriction operator turns a microscopic state into a macroscopic state. The function takes microscopic state (python dictionary) and returns macroscopic state (python dictionary). 

# Installation
Clone the repository to your machine and run 
```
pip install .
```
in the root directory of this project.

# Documentation
I started to work on a documentation using `sphinx`. Currently, the documentation is not published. You can generate a html-file by running
```
make html
```
in the docs folder. (CURRENTLY NOT POSSIBLE)


# Literature
Kevrekidis, I. G., Gear, C. W., & Hummer, G. (2004). Equation‐free: The computer‐aided analysis of complex multiscale systems. AIChE Journal, 50(7), 1346-1355.

Marschler, C., Sieber, J., Berkemer, R., Kawamoto, A., & Starke, J. (2014). Implicit methods for equation-free analysis: convergence results and analysis of emergent waves in microscopic traffic models. SIAM Journal on Applied Dynamical Systems, 13(3), 1202-1238.
