#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 24 18:13:51 2018

@author: paul
"""

import numpy as np
from os import path
from types import NoneType

class state(object):
    def __init__(self,name,purpose):
        assert name is ("macro_state" or "micro_state" or "parameters")
        self.__name = name
        
        assert purpose is ("tmp" or "ref")
        self.__purpose = purpose
        
    @property
    def data(self):
        try:
            return self.__data
        except:
            raise AttributeError("'state' object has no attribute '_state__data'. It has to be set first.")
    
    @data.setter
    def data(self,data_dict):
        assert type (data_dict) is dict
        self.__data = data_dict
        
    def save(self,index=None):
        if self.__purpose == "tmp" and not index == NoneType:
            np.save(path.join(self.__purpose,self.__name+str(index)),self.data)
        else: 
            raise ValueError("keyword argument 'index' must be used if state purpose is 'tmp'")
        
        if self.__purpose == "ref":
            np.save(path.join(self.__purpose,self.__name),self.data)
        
    def load(self,index=None):
        if self.__purpose == "tmp" and not index == NoneType:
            self.data = np.load(path.join(self.__purpose,self.__name+str(index)+".npy")).item()
        else: 
            raise ValueError("keyword argument 'index' must be used if state purpose is 'tmp'")
        
        if self.__purpose == "ref":
            self.data = np.load(path.join(self.__purpose,self.__name+".npy")).item()
        
if __name__ =="__main__":
    data_dict = {"L":np.arange(9,111)}
    
    macro_state = state("macro_state","tmp")
    macro_state.data = data_dict
    macro_state.load(index=1000)