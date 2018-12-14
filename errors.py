#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 11:44:59 2018

@author: paul
"""
import inspect
# =============================================================================
# 
# =============================================================================
class OperatorArgumentError(Exception):
    """Error if the operator has wrong arguments"""
    def __init__(self, operator,operator_args, correct_args):
        msg1 = "Wrong arguments in %s. " %operator.func_name
        msg2 = "%s had the argument names %s. "%(operator.func_name,operator_args)
        msg3 = "%s must have the argument names: %s "%(operator.func_name,correct_args)
        msg = ''.join([msg1,msg2,msg3])
        super(OperatorArgumentError, self).__init__(msg)

def wrong_arguments(operator, correct_args):
    """ Throw an OperatorArgumentError if an operator does not have the correct arguments"""
    operator_args = inspect.getargspec(operator).args
    if operator_args != correct_args:
        raise OperatorArgumentError(operator, operator_args, correct_args) 

# =============================================================================
#         
# =============================================================================
class StateObjectKeyError(Exception):
    """ Error if the provided operator wants to manipulate a variable that was 
    not present in the initial state"""
    def __init__(self,test_key,correct_key_list):
        msg1 = "One operator wanted to assign values to the key %s."%test_key 
        msg2 = "The key was not provided during the initalisation of the stateObject."
        msg = ''.join([msg1,msg2])
        super(StateObjectKeyError, self).__init__(msg)
        
def check_keys(test_key, correct_key_list):
    if test_key not in correct_key_list:
        raise StateObjectKeyError(test_key, correct_key_list)