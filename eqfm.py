# -*- coding: utf-8 -*-
"""
Module for equation-free analysis

TODO: micro and macro state as property
@author: Paul Petersik
"""
import traffic_model as tm
from types import NoneType
import numpy as np
import warnings
import gc
gc.collect()
outpath = "plots/"

class eqfm(object):
    def __init__(self,micro_model,micro_model_parameters,micro_state,macro_state,lifitng_operator,evolution_operator,restriction_operator):
        """
        The class equation free method 
        
        :param points: number of fixed points to be computed
        :param restriction_operator: restriction operator
        :param lifting_operator
        
        :type points: int
        :type points: function
        """
        
        
        self.micro_model = micro_model
        
        assert type (micro_model_parameters) is dict
        self.micro_model_parameters = micro_model_parameters
                
        self.__micro_state = micro_state
        self.__macro_state = macro_state
        
        self.evolution_operator = evolution_operator
        self.restriction_operator = restriction_operator
        self.lifting_operator = lifitng_operator
    
    def get_micro_state(self):
        return self.__micro_state
    
    def get_macro_state(self):
        return self.__macro_state
    
    @staticmethod
    def state(variable_names,dimension):
        variable_dict = {}
        for i in range(len(variable_names)):
            variable_dict[variable_names[i]] = np.zeros(dimension)
        return variable_dict    
    
    
    def lift(self,new_macro_state,new_model_parameters=None):
        """ Lift the macroscopic state into a microscopic state
        """
        self.__micro_state = self.lifting_operator(self,new_macro_state, new_model_parameters)
        self.__macro_state = self.restriction_operator(self,self.__micro_state)
    
    def evolve(self,integration_time):
        """ Evolve the microscopic state
        """
        self.__micro_state = self.evolution_operator(self,integration_time)
    
    def restrict(self):
        """ Restrict the microscopic state. Hence, get the corresponding macroscopic state
        """
        self.__macro_state = self.restriction_operator(self,self.__micro_state)
            
    def compute_reference(self,integration_time):
        """ Compute a reference microscopic state
        """
        self.ref_micro_model_parameters = self.micro_model_parameters
        self.ref_micro_state = self.evolution_operator(self,integration_time,reference=True)
        self.ref_macro_state = self.restriction_operator(self,self.ref_micro_state)

    def compute_macro_time_stepper(self,macro_state_init,tskip,delta,implicit=False):
        """
        macroscopic timestepper
        """
        
        self.lift(macro_state_init)
        self.evolve(tskip)
        self.restrict()
        macro_state_tksip = self.__macro_state
        
        self.evolve(delta)
        self.restrict()
        macro_state_tksip_plus_delta = self.__macro_state
        
        if not implicit:
            explicit_macro_time_stepper = {}
            for key in self.__macro_state.keys():
                explicit_macro_time_stepper[key] = macro_state_tksip_plus_delta[key] - macro_state_tksip[key]  
            return macro_state_tksip, macro_state_tksip_plus_delta, explicit_macro_time_stepper
        
        elif implicit:
            macro_state_init_delta = macro_state_tksip_plus_delta
            macro_state_init_delta2 = {}
            
            for key in self.__macro_state.keys():
                macro_state_init_delta2[key] = 1.1 * macro_state_init_delta[key]
            
            error = {}
            error2 = {}
            derivative_error = {}
            
            for i in range(5):
                self.lift(macro_state_init_delta)
                self.evolve(tskip)
                self.restrict()
                macro_state_tksip = self.__macro_state
                
                for key in self.__macro_state.keys():
                    error[key] =  abs(macro_state_tksip[key] - macro_state_tksip_plus_delta[key])
                    
                self.lift(macro_state_init_delta2)
                self.evolve(tskip)
                self.restrict()
                macro_state_tksip2 = self.__macro_state
                
                for key in self.__macro_state.keys():
                    error2[key] = abs(macro_state_tksip2[key] - macro_state_tksip_plus_delta[key])
                    
                for key in self.__macro_state.keys():
                    derivative_error[key] = (error[key] - error2[key])/(0.1* macro_state_init_delta[key])
                    
                    macro_state_init_delta[key] = macro_state_init_delta[key] - 0.1*error[key]/derivative_error[key]

                    
            implicit_macro_time_stepper = {}
            for key in self.__macro_state.keys():
                implicit_macro_time_stepper[key] = (macro_state_init_delta[key] -  macro_state_init[key])
            
            return macro_state_init, macro_state_init_delta, implicit_macro_time_stepper
    
    def compute_macro_time_derivative(self,macro_time_stepper,delta):
        """
        time derivative of macroscopic state evolution
        """
        F = {}
        for key in self.__macro_state.keys():
            F[key] = macro_time_stepper[key]/delta
        return F
        
    def project_macro_state(self,macro_state_init,tskip,delta,Delta_t,implicit=False):
        """
        projective integration
        """
        macro_state_tksip, macro_state_tksip_plus_delta, macro_time_stepper = self.compute_macro_time_stepper(macro_state_init,tskip,delta,implicit)
        F = self.compute_macro_time_derivative(macro_time_stepper,delta)
        
        projected_macro_state = {}
        if Delta_t>0:
            for key in self.__macro_state.keys():
                projected_macro_state[key] = macro_state_tksip_plus_delta[key] + Delta_t * F[key]
        
        if Delta_t<0:
            for key in self.__macro_state.keys():
                projected_macro_state[key] = macro_state_tksip[key] + Delta_t * F[key]
        
        return projected_macro_state
            
    def projective_integration(self,tskip,delta,Delta_t,iterations,macro_state_init = None,implicit=False,verbose=False):
        """ coarse time stepper using extrapolation of the macropscopic state
        from the microscopic model"""
        
        self.tskip = tskip
        self.delta = delta
        if type(macro_state_init) != NoneType:
            self.__macro_state  = macro_state_init
        
        if Delta_t<0 and (tskip + Delta_t)>0:
            warnings.warn("backwards integration might be ineffective because tksip>Delta_t",Warning)
        
        for i in range(iterations):
            self.__macro_state = self.project_macro_state(self.__macro_state,tskip,delta,Delta_t,implicit)
            if verbose:
                print self.__macro_state
    



# =============================================================================
# =============================================================================
# # Example traffic model
# =============================================================================
# =============================================================================

if __name__ == "__main__":
        
    # initalize micro model
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
    micro_state = eqfm.state(micro_var_names,number_of_cars)
    
    # initialize macro state
    macro_var_name = ["standard_deviation_headway","standard_deviation_velocity"]
    macro_dim = 1
    macro_state = eqfm.state(macro_var_name,macro_dim)  
    macro_state["standard_deviation_headway"] = 5.
    macro_state["standard_deviation_velocity"] = 4.
    
    # define lifting operator
    def lifting_operator(self,new_macro_state,new_micro_model_parameters=None):
        
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
        
        micro_model = self.micro_model(self.micro_model_parameters)
        
        if not reference:
            micro_state = self.get_micro_state()
            x_init = micro_state["position"]
            dot_x_init = micro_state["velocity"]
            ddot_x_init = micro_state["acceleration"]
            micro_model.initCars(x_init=x_init,dot_x_init=dot_x_init,ddot_x_init=ddot_x_init)
        
        else:
            micro_model.initCars()
            micro_state = {}
        
        micro_model.integrate()
        
        micro_state["position"] = micro_model.x[:,-1]
        micro_state["velocity"] = micro_model.dot_x[:,-1]
        micro_state["acceleration"] = micro_model.ddot_x[:,-1]
        micro_state["headway"] = micro_model.Delta_x[:,-1]
        
        del(micro_model)
        gc.collect()
        
        return micro_state
      
    # define restriction operator
    def restriction_operator(self,micro_state):
        macro_state = {}
        macro_state["standard_deviation_headway"] = np.std(micro_state["headway"])
        macro_state["standard_deviation_velocity"] = np.std(micro_state["velocity"])
        #macro_state["L"] = self.micro_model_parameters["L"]
        return macro_state
    
    
    
    # initialize the equation free model

    model = eqfm(traffic_model,
                      traffic_model_parameters,
                      micro_state,
                      macro_state,
                      lifting_operator,
                      evolution_operator,
                      restriction_operator)
    
# =============================================================================
# =============================================================================
# # Run application
# =============================================================================
# =============================================================================
    
    print "============COMPUTE REFERENCE =============="
    model.compute_reference(150)
    print model.ref_macro_state
    
    print "========PROJECTIVE INTEGRATION========"
    model.projective_integration(5,20,30,1,implicit=True,verbose=True)
    
    del(model)


gc.collect()