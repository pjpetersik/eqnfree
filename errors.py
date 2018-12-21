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

def wrong_arguments(operator, expected_args):
    """ Throws an OperatorArgumentError if an operator does not have the correct arguments/argument names 
    needed for the module to be compatible with the from the user provided operators
    
    :param operator: the operator function that shell be tested for the expected arguments
    :param expected_args: a list of strings of the expected arguments of the operator function
    
    :type operator: function
    :type expected_args: list of strings
    """
    operator_args = inspect.getargspec(operator).args
    if operator_args != expected_args:
        raise OperatorArgumentError(operator, operator_args, expected_args) 

# =============================================================================
#         
# =============================================================================
class StateObjectKeyError(Exception):
    """ Error if the provided operator wants to manipulate a variable that was 
    not present in the initial state"""
    def __init__(self,test_key,correct_key_list, category, purpose):
        msg1 = "Used key: '%s' in stateObject (category: %s, purpose: %s)" %(test_key,category,purpose)
        msg2 = " is not in list of initialised  keys: %s" %correct_key_list
        msg = ''.join([msg1,msg2])
        super(StateObjectKeyError, self).__init__(msg)
        
def check_keys(test_key_list, expected_key_list, category, purpose):
    """ Throws an error if the key used for a state dictionary is not in the list of 
    keys of the initialised stateObject (expected_key_list). Note, ref and tmp stateObject 
    of the same category ("micro","macro","parameters") have the same 'expected_key_list'.
    
    :param test_key_list: list of keys that should be tested against the expected_key_list
    :param expected_key_list: the expected_key_list
    
    :type test_key_list: list 
    :type expected_key_list: list 
    """
    
    for test_key in test_key_list:
        if test_key not in expected_key_list:
            raise StateObjectKeyError(test_key, expected_key_list, category, purpose)