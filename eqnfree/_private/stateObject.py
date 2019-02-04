#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 27 18:31:27 2018

@author: paul
"""
from types import NoneType
import numpy as np
import os
from errors import  check_keys, check_for_dictType, check_purpose, check_category

class stateObject(object):
    """
    Holds the values for the micro, macro or parameter state of a equation-free
    model. The state values can be hold for 2 purposes either for a reference case ("ref") 
    or for a temporary case ("tmp"). Furthermore, the category of a state can be 
    either "micro","macro" or"parameters" refering to microscopic, macroscopic and parameter state.
    
    """
    def __init__(self,category,purpose,keys, data=None):
        """
        :type category: str
        :param category: Either "micro","macro" or"parameters".
        
        :type purpose: str
        :param purpose: Either "tmp" or "ref".
        
        :type data: dict
        :param data: The values of the considered state.
        
        """
        
        check_category(category)
        check_purpose(purpose)
        check_for_dictType(data, category, purpose)
        
        self.__category = category
        self.__purpose = purpose
        self.__keys = keys
        
        if type(data) is dict:
            self.__data = data
        elif type(data) is NoneType:
            self.__data = {}
    
    @property
    def variableDict(self):
        return self.__data
    
    @variableDict.setter
    def variableDict(self,new_data):

        check_for_dictType(new_data,self.__category,self.__purpose)
        
        test_key_list = new_data.keys()
        check_keys(test_key_list, self.__keys, self.__category,self.__purpose)
        
        self.__data = new_data
    
    @property
    def category(self):
        return self.__category
    
    @property
    def purpose(self):
        return self.__purpose
    
    def save(self,index=None):
        """
        Saves the data dictionary to a NPY-file to a distinct location.
        
        :type index: int
        :param index: The index of the current time step for which the data should be  saved.
        
        """
        outputFolder = self.purpose
        if self.purpose =="tmp":
            assert index!=NoneType          
            outputFile = self.category + str(index)+".npy"
           
        elif self.purpose =="ref":
            outputFile = self.category +".npy"  
        outputPath = os.path.join(outputFolder,outputFile)
        
        if not os.path.exists(outputFolder):
            os.mkdir(outputFolder)
        
        np.save(outputPath,self.__data)
    
    def load(self,index=None):
        """
        Loads data from previous runs that where saved using the method stateObject.save()
        
        :type index: int
        :param index: If the category of a stateObject is "tmp" an index has to be provided for which the data should be loaded.
        
        """
        inputFolder = self.purpose
        if self.purpose =="tmp":
            assert index!=NoneType          
            inputFile = self.category + str(index)+".npy"
           
        elif self.purpose =="ref":
            inputFile = self.category +".npy"
        
        inputPath = os.path.join(inputFolder,inputFile)
        
        inputDict = np.load(inputPath).item()
        for key in inputDict.keys():
            self.__data[key] = inputDict[key]