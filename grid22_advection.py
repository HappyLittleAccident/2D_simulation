# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 20:47:52 2022

@author: Marek
"""

import numpy as npy
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

plt.close('all')

def system(t,y,a,b,d,l,w,gp,gn):
    np = y[:4]
    nn = y[4:]
    
    dot_np = [
        #dot_np[0]
        -a*((np[2]-nn[2])*(-(np[1]-np[0])/l)+(np[1]-nn[1])*((np[2]-np[0])/w) #outgoing
              -(np[3]-nn[3])*((np[0]-np[1])/l)-(np[3]-nn[3])*(-(np[0]-np[2])/w)) #ingoing from other squares
        + b*(np[0]+nn[0]) - d*nn[0]*np[0] + gp[0],
        #dot_np[1]        
        -a*((np[3]-nn[3])*((np[0]-np[1])/l)+(np[0]-nn[0])*(-(np[3]-np[1])/w)
              -(np[2]-nn[2])*(-(np[1]-np[0])/l)-(np[2]-nn[2])*((np[1]-np[3])/w))        
        + b*(np[1]+nn[1]) - d*nn[1]*np[1] + gp[1],
        #dot_np[2]        
        -a*((np[0]-nn[0])*((np[3]-np[2])/l)+(np[3]-nn[3])*(-(np[0]-np[2])/w)
              -(np[1]-nn[1])*(-(np[2]-np[3])/l)-(np[1]-nn[1])*((np[2]-np[0])/w))        
        + b*(np[2]+nn[2]) - d*nn[2]*np[2] + gp[2],
        #dot_np[3]      
        -a*((np[1]-nn[1])*(-(np[2]-np[3])/l)+(np[2]-nn[2])*((np[1]-np[3])/w)
              -(np[0]-nn[0])*((np[3]-np[2])/l)-(np[0]-nn[0])*(-(np[3]-np[1])/w))        
        + b*(np[3]+nn[3]) - d*nn[3]*np[3] + gp[3]
        
        ]
    dot_nn = [        
        #dot_nn[0]
        -a*((np[2]-nn[2])*(-(nn[1]-nn[0])/l)+(np[1]-nn[1])*((nn[2]-nn[0])/w)
              -(np[3]-nn[3])*((nn[0]-nn[1])/l)-(np[3]-nn[3])*(-(nn[0]-nn[2])/w)) #ingoing from other squares
        + b*(np[0]+nn[0]) - d*nn[0]*np[0] + gn[0],
        #dot_nn[1]        
        -a*((np[3]-nn[3])*((nn[0]-nn[1])/l)+(np[0]-nn[0])*(-(nn[3]-nn[1])/w)
              -(np[2]-nn[2])*(-(nn[1]-nn[0])/l)-(np[2]-nn[2])*((nn[1]-nn[3])/w)) 
        + b*(np[1]+nn[1]) - d*nn[1]*np[1] + gn[1],
        #dot_nn[2]        
        -a*((np[0]-nn[0])*((nn[3]-nn[2])/l)+(np[3]-nn[3])*(-(nn[0]-nn[2])/w)
              -(np[1]-nn[1])*(-(nn[2]-nn[3])/l)-(np[1]-nn[1])*((nn[2]-nn[0])/w))
        + b*(np[2]+nn[2]) - d*nn[2]*np[2] + gn[2],
        #dot_nn[3]      
        -a*((np[1]-nn[1])*(-(nn[2]-nn[3])/l)+(np[2]-nn[2])*((nn[1]-nn[3])/w)
              -(np[0]-nn[0])*((nn[3]-nn[2])/l)-(np[0]-nn[0])*(-(nn[3]-nn[1])/w)) 
        + b*(np[3]+nn[3]) - d*nn[3]*np[3] + gn[3]

        ]

    return npy.ravel([dot_np,dot_nn])

gp = [1,0,0,0]
gn = [0,0,1,0]
w=1
l=1
a=2
b=1
d=4     
        #np     #nn
       #1 2 3 4 1 2 3 4
y0s = [[0,0,0,0,0,0,0,0],
       [0,1,1,0,1,0,0,1],
       [1,0,0,0,0,1,0,0]
       ]
for y0 in y0s:   
    
    sol = solve_ivp(system,(0,1e4),y0,args=(a,b,d,l,w,gp,gn))
    
    np_final = sol.y[:4,-1].reshape(2,2)
    nn_final = sol.y[4:,-1].reshape(2,2)
    circ = (sol.y[:4,-1].reshape(2,2)-sol.y[4:,-1].reshape(2,2))

    scale = npy.max(circ)
    
    fig_22, ax_22 = plt.subplots(1, 1)
    im = ax_22.imshow(circ, aspect='auto', cmap='seismic', vmin=-scale, vmax=scale)
    fig_22.colorbar(im)
    fig_22.tight_layout()






