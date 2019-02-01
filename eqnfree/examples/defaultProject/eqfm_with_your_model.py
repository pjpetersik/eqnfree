#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 25 17:28:46 2018

@author: paul
"""

import numpy as np
import gc
from eqfm import eqfModel
# =============================================================================
# import your model
# =============================================================================
# This model should contain methods for the evolution/integration of the microscopic
# model e.g. your_model.integrate() that you are going to use later
from your_model_modlue import your_model

# =============================================================================
# # Put all parameters of your model in to a dictionary
# =============================================================================
parameters = { "parameter1_name": 1 }

# =============================================================================
# # initialize micro state
# =============================================================================
# names of micro states
micro_var_names = ["micro_state1_name","micro_state2_name"]
# dimension of the micro states
micro_dim = 60
# helper function to initialize a dictionary with keys
# from the list "micro_var_names" and numpy zero arrays with dimension "micro_dim"
initial_micro_state = eqfModel.state(micro_var_names,micro_dim)

# =============================================================================
# # initialize macro state
# =============================================================================
# Note that currently the methods of the eqfModel module are just applied to one 
# macro state that has dimension 1!

# names of macro states
macro_var_name = ["macro_state1_name","macro_state2_name"]
# dimension of the macro states
macro_dim = 1
# helper function to initialize a dictionary
initial_macro_state = eqfModel.state(macro_var_name,macro_dim)  

# =============================================================================
# =============================================================================
# # Define operators
# =============================================================================
# =============================================================================

# The Lifting, evolution and restriction operators will be defined in the following.
# They will be used in the generation of the eqfModel instance. And will be 
# available as methods of the eqfModel instance. This is way, they have the self
# argument. Hence, understand the operators as methods of the eqfModel instance.
# With this you can use the reference states in the lifting operation and the the 
# current micro states in the restriction operation.

# =============================================================================
# # define lifting operator
# =============================================================================
# arguments MUST stay as they are in here
def lifting_operator(self,new_macro_state,new_micro_model_parameters=None):
    # make some checks (if wanted)
    assert new_macro_state["macro_state1_name"]>0
    
    # change parameter value/the state the model should be lifted into    
    if new_micro_model_parameters is not None:
        self.micro_model_parameters["parameter1_name"] = new_micro_model_parameters["parameter1_name"]    
    
    # use the reference micro state that is saved in the eqfModel instance
    micro_state1_ref = self.ref_micro_state["micro_state1_name"]
    
    # use the referencemacro state that is saved in the eqfModel instance
    macro_state1_ref = self.ref_micro_state["macro_state1_name"]
    
    # use the referencemacro state that is saved in the eqfModel instance
    parameter1_ref = self.ref_micro_model_parameters["parameter1_name"]
    
    # Now mace add your lifting operator
    micro_state1 = 1*2*3
    # ....
    
    # save the micro_state values to dictionary
    micro_state = {}
    micro_state["micro_state1_name"] = micro_state1
    
    # return the dictionary of the lifted state!
    return micro_state

# =============================================================================
# # define evolution operator
# =============================================================================
#  AGAIN: Dont change the arguments of the function
def evolution_operator(self,integration_time,reference = False):
    # change the some parameters of the model.
    self.micro_model_parameters["parameter1_name"] = 2
    
    # It might be useful to have the integration time as value in the parameter dictionary (but not  needed)
    self.micro_model_parameters["tmax"] = integration_time
    
    # make some changes if the you are NOT computing a reference state
    if not reference:
        #....some code
        micro_state_dict = self.micro_state
    
    # make some changes if the you ARE computing a reference state
    else:
        # some initalisation code 
        self.micro_model.init() # not necessary needed to be done with the init() method
        micro_state_dict = {}
    
    # compute the microscopic evolution/integration
    #.....some code such as
    micro_state1 = self.micro_model.integrate() # the evolution does not need to be done with a mirco_model method (but can be)
    
    micro_state_dict["micro_state1_name"] = micro_state1

    # clean up stuff to free memory
    gc.collect()
    
    return micro_state_dict

# =============================================================================
# # define restriction operator
# =============================================================================
#AGAIN: dont  change the arguemnts
def restriction_operator(self,micro_state):
    macro_state_dict = {}
    
    # compute the macro_state from the micro state
    #... e.g. the mean and std:
    macro_state_dict["macro_state1_name"] = np.mean(micro_state["micro_state1_name"])
    macro_state_dict["macro_state2_name"] = np.std(micro_state["micro_state2_name"])
    
    # return  macro state dictionary
    return macro_state_dict

# =============================================================================
# =============================================================================
# # # Generate a eqfModel instance
# =============================================================================
# =============================================================================
model = eqfModel(your_model,
                  parameters,
                  initial_micro_state,
                  initial_macro_state)

# set the lift, evolve, restric operators
model.setEqfmOperators(lifting_operator,
                  evolution_operator,
                  restriction_operator)

# set tksip, delta and boolean for implicit or explicit method
tskip = 10
delta = 200
implicit = True # implicit method-> True, explicit method -> False

model.setEqfmParameters(tskip,delta,implicit)

# =============================================================================
# =============================================================================
# # Run bifurcation analysis
# =============================================================================
# =============================================================================
number_fixed_points = 400
# there are a lot of keyword arguments that you can use and most likely need to use to make
# the code working for you (check the doc string)
model.bifurcation_analysis("parameter1_name","macro_state1_name", number_fixed_points)

# =============================================================================
# =============================================================================
# # Run projective integration
# =============================================================================
# =============================================================================
delta_t = 35.
number_of_iterations = 1000
model.projective_integration(delta_t,number_of_iterations,"macro_state1_name")

gc.collect()