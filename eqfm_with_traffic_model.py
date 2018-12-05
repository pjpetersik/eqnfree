#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 25 17:28:46 2018

@author: paul
"""
import traffic_model as tm
import numpy as np
import gc
from eqfm import eqfModel
from types import NoneType
    
# micro model class
traffic_model = tm.model.from_dictionary

# initalize model parameters
traffic_model_parameters = {}
traffic_model_parameters["N"] = 22
traffic_model_parameters["L"] = 250.
traffic_model_parameters["dt"] = 1./3.
traffic_model_parameters["tmax"] = 25.
traffic_model_parameters["xpert"] = 1.*np.sin(2*np.pi/float(traffic_model_parameters["N"])
                                    * np.arange(traffic_model_parameters["N"]))
# initialize micro state
micro_var_names = ["position","velocity","acceleration","headway"]
number_of_cars = 22
micro_state = eqfModel.state(micro_var_names,number_of_cars)

# initialize macro state
macro_var_name = ["standard_deviation_headway","standard_deviation_velocity"]
macro_dim = 1
macro_state = eqfModel.state(macro_var_name,macro_dim)  
macro_state["standard_deviation_headway"] = 3.
macro_state["standard_deviation_velocity"] = 5.

# define lifting operator
def lifting_operator(self,new_macro_state,new_micro_model_parameters=None):
    assert new_macro_state["standard_deviation_headway"]>0
    assert new_macro_state["standard_deviation_velocity"]>0
    
    std_Dx = new_macro_state["standard_deviation_headway"]
    std_dotx = new_macro_state["standard_deviation_velocity"]
    
    if type(new_micro_model_parameters) != NoneType:
        L = new_micro_model_parameters["L"]
        self.micro_model_parameters["L"] = new_micro_model_parameters["L"]    
    else:
        L = self.micro_model_parameters["L"]
        
    Dx_ref = self.ref_micro_state["headway"]
    dotx_ref = self.ref_micro_state["velocity"]
    
    std_Dx_ref = self.ref_macro_state["standard_deviation_headway"]
    std_dotx_ref = self.ref_macro_state["standard_deviation_velocity"]
    
    L_ref = self.ref_micro_model_parameters["L"]
    
    x = np.zeros_like(self.ref_micro_state["position"])
    dotx = np.zeros_like(x)
    ddotx = np.zeros_like(x)
    Dx = np.zeros_like(x)
    
    Dx = std_Dx/std_Dx_ref * (Dx_ref - Dx_ref.mean()) + L/L_ref * Dx_ref.mean()
        
    x[0] = 0
    x[1:] = np.cumsum(Dx[:])[:-1]
        
    dotx[:] =  std_dotx/std_dotx_ref * (dotx_ref - dotx_ref.mean()) + dotx_ref.mean()
    ddotx[:] = 0
    
    micro_state = {}
    micro_state["position"] = x
    micro_state["velocity"] = dotx
    micro_state["acceleration"] = ddotx
    micro_state["headway"] = Dx
    
    return micro_state

# define evolution operator
def evolution_operator(self,integration_time,reference = False):
    self.micro_model_parameters["tmax"] = integration_time
    self.micro_model.update_parameters(self.micro_model_parameters)
    
    if not reference:
        micro_state = self.micro_state
        x_init = micro_state["position"]
        dot_x_init = micro_state["velocity"]
        ddot_x_init = micro_state["acceleration"]
        self.micro_model.initCars(x_init=x_init,dot_x_init=dot_x_init,ddot_x_init=ddot_x_init)
    
    else:
        self.micro_model.initCars()
        micro_state = {}
    
    self.micro_model.integrate()
    
    micro_state["position"] = self.micro_model.x[:,-1]
    micro_state["velocity"] = self.micro_model.dot_x[:,-1]
    micro_state["acceleration"] = self.micro_model.ddot_x[:,-1]
    micro_state["headway"] = self.micro_model.Delta_x[:,-1]
    
    #del(micro_model)
    gc.collect()
    return micro_state
  
# define restriction operator
def restriction_operator(self,micro_state):
    macro_state = {}
    macro_state["standard_deviation_headway"] = np.std(micro_state["headway"])
    macro_state["standard_deviation_velocity"] = np.std(micro_state["velocity"])
    return macro_state

# initialize the equation free model

model = eqfModel(traffic_model,
                  traffic_model_parameters,
                  micro_state,
                  macro_state)

model.setEqfmOperators(lifting_operator,
                  evolution_operator,
                  restriction_operator)

model.setEqfmParameters(10,20,False)
# =============================================================================
# =============================================================================
# # Run application
# =============================================================================
# =============================================================================


model.bifurcation_analysis("L","standard_deviation_headway",100,dmacro = 0.1,s=[0.1,10],rerun=False)

model.projective_integration(35.,100,"standard_deviation_headway")
#%%
import matplotlib.pyplot as plt
plt.scatter(model.fixed_points["L"],model.fixed_points["standard_deviation_headway"])
#del(model)
gc.collect()