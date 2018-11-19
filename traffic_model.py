#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 11:43:14 2018

@author: paul
"""
import numpy as np
from tensorflow import keras as kr
from types import NoneType 
import gc

class model(object):
    def __init__(self,N=None,L=None,tmax=None,dt=None,xpert=0):         
        self.N =  N
        self.L =  L
        self.distance = np.arange(0,self.L,1) #array for distance
        self.xpert = xpert
        
        self.tmax  = tmax
        self.dt    = dt
        self.iters = abs(int(self.tmax/self.dt))
        self.time  = np.arange(0,self.tmax,self.dt)

        self.loaded_tf_model = kr.models.load_model('model/model.h5')
        self.meanDx,self.stdDx,self.meandotx,self.stddotx,self.n_leading_cars,self.rnn = np.loadtxt('model/model_parameter.txt')
        self.n_leading_cars  = int(self.n_leading_cars)
        
    def update_parameters(self,parameters):         
        self.N =  parameters["N"]
        self.L =  parameters["L"]
        self.xpert = parameters["xpert"]
        self.tmax  = parameters["tmax"]
        self.dt    = parameters["dt"]
        
        self.distance = np.arange(0,self.L,1) #array for distance
        self.iters = abs(int(self.tmax/self.dt))
        self.time  = np.arange(0,self.tmax,self.dt)
    
    def __del__(self):
        """
        free memory
        """
        kr.backend.clear_session()
                
    @classmethod
    def from_dictionary(cls,parameters):
        return cls(N=parameters["N"], L=parameters["L"], tmax=parameters["tmax"],dt=parameters["dt"],xpert=parameters["xpert"])
        
    def initCars(self, x_init=None, dot_x_init = None, ddot_x_init=None):
        """
        initialise 0th time step
        """  
        
        self.x       = np.zeros(shape=(self.N,self.iters)) # position
        self.dot_x   = np.zeros(shape=(self.N,self.iters)) # velocity
        self.ddot_x  = np.zeros(shape=(self.N,self.iters)) # acceleration
        self.Delta_x = np.zeros(shape=(self.N,self.iters)) # headway
        
        if type(x_init)==NoneType:
            self.x[:,0] = np.linspace(0,self.L,self.N) 
            self.x[:,0] = self.x[:,0] + self.xpert
        else:
            self.x[:,0] = x_init
            
        if type(dot_x_init)==NoneType:
            self.dot_x[:,0]  = self.meandotx
        else:
            self.dot_x[:,0]  = dot_x_init
            
        if type(ddot_x_init)==NoneType:
            self.ddot_x[:,0] = 0.
        else:
            self.ddot_x[:,0] = ddot_x_init
            
        self.Delta_x[:,0]   = self.headway(self.x[:,0],self.L)
        
        
            
        
    def integrate(self):
        """
        Integrate the model using a fortran or a python kernel 
        """
        for i in range(0,self.iters-1):
            self.integration_procedure(i)
            
    def integration_procedure(self,i):
        """
        RK4 integration scheme
        """
        h = self.dt
        k1 = self.acceleration(self.Delta_x[:,i],self.dot_x[:,i],self.Delta_x[:,i-1],self.dot_x[:,i-1])
        self.dot_x[:,i+1] = self.dot_x[:,i] + k1*h/2
        
        k2 = self.acceleration(self.Delta_x[:,i],self.dot_x[:,i+1],self.Delta_x[:,i-1],self.dot_x[:,i-1])
        
        self.dot_x[:,i+1] = self.dot_x[:,i] + k2*h/2
        k3 = self.acceleration(self.Delta_x[:,i],self.dot_x[:,i+1],self.Delta_x[:,i-1],self.dot_x[:,i-1])
        
        self.dot_x[:,i+1] = self.dot_x[:,i] + k3*h
        k4 = self.acceleration(self.Delta_x[:,i],self.dot_x[:,i+1],self.Delta_x[:,i-1],self.dot_x[:,i-1])
        
        self.ddot_x[:,i+1] = k1 
        
        self.dot_x[:,i+1] = self.dot_x[:,i] + h/6. * (k1 + 2*k2 + 2*k3 + k4) 
        
        # just allow postive velocities
        self.dot_x[self.dot_x[:,i+1]<0.,i+1] = 0.
        
        self.x[:,i+1]      = self.x[:,i] + self.dot_x[:,i+1] * h
        
        self.x[:,i+1]      = self.x[:,i+1]%self.L

        # Diagnostics
        self.Delta_x[:,i+1]   = self.headway(self.x[:,i+1],self.L)

    def acceleration(self,Delta_x,dot_x,prev_Delta_x,prev_dot_x):
        """
        returns the acceleration of cars based on an ANN
        """
        def shape(Delta_x,dot_x):
            Dx = Delta_x.reshape((len(Delta_x),1))
            dotx = dot_x.reshape((len(Delta_x),1))
                    
            Dx = (Dx - self.meanDx)/self.stdDx
            dotx = (dotx - self.meandotx)/self.stddotx
            
            X = np.concatenate((Dx[0:self.n_leading_cars,:].T,dotx[0:self.n_leading_cars,:].T),axis=1)
    
            X_all = np.zeros((self.N,len(X[0])))
            
            # for car number 0
            X_all[0,:] = X
            
            for i in range(self.N-1):    
                Dx = np.roll(Dx,-1,axis=0)
                dotx = np.roll(dotx,-1,axis=0)            
                X_all[i+1,:] =  np.concatenate((Dx[0:self.n_leading_cars,:].T,dotx[0:self.n_leading_cars,:].T),axis=1)
            
            if self.rnn:
                X_all = X_all.reshape(X_all.shape[0],1,X_all.shape[1])
            return X_all
        
        if self.rnn:
            X = shape(Delta_x,dot_x)
#            X_prev = shape(prev_Delta_x,prev_dot_x)
#            
#            X = np.concatenate((X,X_prev),axis=1)
        else:
            X = shape(Delta_x,dot_x)
            
        # predict acceleration
        ddotx = self.loaded_tf_model.predict(X)
        
        return ddotx[:,0]
    
    def headway(self,x,L):
        Dx = np.zeros(self.N)
        Dx[:-1] = ((x[1:] - x[:-1])+L)%L
        Dx[-1] = ((x[0] - x[-1])+L)%L
        return Dx 
    
    