# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 14:56:40 2020

@author: Admin
"""

import numpy as np
import math as m
import matplotlib.pyplot as plt


#argument as follows pixels, exposure time, quantum efficiency, area of pixel in micron, photon flux 
class noise():
    NT=0
    SNR=0
    s=0
    nd=0
    nw=0
    nr=0
    nf=0
    sk=0
    shot=0
    nd2=0
    nrr = 0
    def get(*args):
        
        #pixels
        pixel=args[0]
        #exposure time
        t=args[1]
        #quantum efficiency
        qe = args[2]/100
        
        #dark noise
        
        #dark current  nanoamperes/cm2 at T = 300K  0-3
        I = 0.018
        
        #area of a pixel in cm2
        # A= 6.25 * np.power(10,-6.0)
        A = args[3]*np.power(10,-8.0)
        
        #dark current at T = 300 K
        #I = n*1.6*np.power(10,-19.0)/A
        
        #temp. in K
        T = 193
       
        
        #boltzmann constant in eV K-1
        k = 8.617333262145*np.power(10,-5.0)
        
        #eV
        E = 1.1557 - ((7.021*np.power(10,-4.0)*np.power(T,2.0))/(1108+T))
        
        #e/pixel/second
        ND = 2.5*np.power(10,15.0)*A*I*np.power(T,1.5)*m.exp(-E/(2*k*T))
       
        noise.nd = ND*t*pixel
        

        
        
        
        
        #shot noise  photons/pixel/second + Q.E.
        #signal or photon flux
        S = args[4]
        
        NS = np.power(S,0.5)
        # print(NS)
        #Quantum efficiency
        noise.shot = S*t*qe
        
        sensitivity = 1*np.power(10,-6.0)
        #readout noise rms eletrons/pixel
        #white noise R = 200 o 20K
        B=100000000
        R=2000
        k=1.380649*np.power(10,-23.0)
        NW = np.power(4*k*T*B*R,0.5)
        
        #in electrons
        noise.nw = NW/sensitivity
        
        
        #reset noise
        B=100000000
        R=500
        k=1.380649*np.power(10,-23.0)
        NR = np.power(4*k*T*B*R,0.5)
    
        #in electrons
        noise.nr = NR/sensitivity
        
        
        #flicker noise approx. nV/rtHz
        NF = 100*np.power(10,-9.0)
        noise.nf = NF/sensitivity
        
        #total readout noise
        nrr = noise.nr+noise.nf+noise.nw
        noise.nrr  = nrr*nrr*pixel
        
        
        
        #sky back noise
        mv=21
        noise.sk  = 180*180*3.14*np.power(10,-0.4*mv)*8.8*np.power(10,9.0)*qe*t*1.2*1.2*3.14*pixel/(1024*1024)
        
        
        #total noise
        noise.NT = np.power(noise.nrr +noise.shot + noise.nd + noise.sk,0.5)

        # print(noise.NT)
        
        #signal to noise ratio
        noise.s = noise.shot
        
        noise.SNR = noise.s/noise.NT
        # print(noise.SNR)
        











#dark noise
# import itertools
# import seaborn as sns

# palette = itertools.cycle(sns.color_palette())
# plt.figure(figsize=(16,16))      
# darklist =[]    
# time = [] 
# T=193
# for t in range(0,1000,100):
   
#     noise.get(1,t,60,169,10000,T)
#     darklist.append(noise.nd2)
#     time.append(t)
    
#     # plt.legend()
#     plt.plot(time,darklist,'r-',color=next(palette),label=str(T) )
#     plt.xlabel("Exposure time (s)",fontsize = 20)
#     plt.xticks(fontsize=18)
#     plt.ylabel("Dark noise (e- per pixel)",fontsize =20)
#     plt.yticks(fontsize=18)
        
# # plt.legend(title='Temperature'+ 'K')    
# plt.show()
# print(darklist)
#argument as follows pixels, exposure time, quantum efficiency, area of pixel in micron, photon flux      

import itertools
import seaborn as sns

palette = itertools.cycle(sns.color_palette())
plt.figure(figsize=(16,16))      
snrl =[]    
time = [] 
s = 450000000

for t in range(0,10000,100):
   
    noise.get(s,t,60,169,20000)
    snrl.append(noise.SNR)
    time.append(t)
    
    
    plt.plot(time,snrl,'r-')
    plt.xlabel("Exposure time (s)",fontsize = 20)
    plt.xticks(fontsize=18)
    plt.ylabel("SNR",fontsize =20)
    plt.yticks(fontsize=18)
    
plt.legend(title='7500 Photon influx (e/sec/pixel)')  
plt.show()