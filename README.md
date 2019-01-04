# Equation-free modeling
In this repository I work on developing a module that contains code from the equation-free modeling developed by Kevrekidis et al. (2004) and refined by Marschler et al. (2014). Up until now, my focus was on developing a method for a bifurcation analysis using the equation free method (as described in Marschler et al. (2014)). In the example directory examples for the application of the equation-free modeling are given.

NOTE: This is work in progress. But the current version is working for the examples given in the examples folder. I am very happy over any feedback and contribution to the project!

# General concept
Writing the code for the equation-free methods (e.g. bifurcation analysis in in Marschler et al. (2014)) can be a tedious work. In this project, I want to develop a framework such that equation-free methods are easily applicable to any microscopic model (e.g. traffic model, puplic goods games, etc). The eventual aim is, that a user just needs to provide the lifting, evolution and restriction operator and the parameter settings to make use of methods such as bifurcation analysis or projective integration. 

# Literature
Kevrekidis, I. G., Gear, C. W., & Hummer, G. (2004). Equation‐free: The computer‐aided analysis of complex multiscale systems. AIChE Journal, 50(7), 1346-1355.

Marschler, C., Sieber, J., Berkemer, R., Kawamoto, A., & Starke, J. (2014). Implicit methods for equation-free analysis: convergence results and analysis of emergent waves in microscopic traffic models. SIAM Journal on Applied Dynamical Systems, 13(3), 1202-1238.

