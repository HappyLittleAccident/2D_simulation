# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 14:33:45 2022

@author: Marek

osageruv plyn !!!!
"""

import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
import time

plt.close('all')

def vortices(t,y,a,b,d,g,gs):
    n,s = y
    return np.array([(a+b)*n  - 0.5*d*n**2*(1-s**2) + g ,
                     -2*b*s + 0.5*d*s*n*(1-s**2) + (g/n)*(-s+gs)])

def populations(t,y,a,b,d,g,gs):
    n1,n2 = y
    return np.array([a*n1 + b*n2 - d*n1*n2 + g , a*n2 +b*n1 - d*n1*n2 + gs])







n02 = np.linspace(0.01,1,15)
n01 = np.linspace(0.01,1,15)

s0s = np.linspace(-1,1,10)
n0s = np.linspace(0.1,3,15)

def stability(t,y,a,b,d,g,gs):
    n,s = y
    if n<-4:
        return 0
    else:
        return np.abs((a+b)*n - 0.5*d*n**2*(1-s**2) + g) + np.abs(-2*b*s+ 0.5*d*s*n*(1-s**2) + (g/n)*(gs-s))-4e-3

def populations_stability(t,y,a,b,d,g,gs):
    n1,n2 = y
    return np.sum(np.abs(np.array([a*n1 -b*n2 - d*n1*n2 + g , a*n2 +b*n1 - d*n1*n2 + gs])))-4e10

stability.terminal = True

a = -1
b = 1
d = 4
g = 1
gs = 1
stable=[]


plt.figure()
for s0 in s0s:
    for n0 in n0s:
        t1 = time.time()
        result = integrate.solve_ivp(fun=vortices,t_span=(0,1e1),y0=np.array([n0,s0]),args = (a,b,d,g,gs),events = stability)
        t2 = time.time()
        print('Integration took: {:.3f}s'.format(t2-t1))
        if result.y[1][-1]>0:
            color = 'orange'
        else:
            color = 'blue'
        plt.plot(result.y[0],result.y[1],c=color)
        stable.append([result.y[0][-1],result.y[1][-1]])
        
stable = np.array(stable)

plt.plot(stable[-1,0],stable[-1,1],'k*',ms=15)
plt.plot(stable[0,0],stable[0,1],'k*',ms=15)
plt.xlabel('n')
plt.ylabel('s')

plt.xlim(0,3)

# =============================================================================
# plt.figure()
# for s0 in s0s:
#     for n0 in n0s:
#         t1 = time.time()
#         result = integrate.solve_ivp(fun=populations,t_span=(0,1e1),y0=np.array([n0,s0]),args = (a,b,d,g,gs),events = populations_stability)
#         t2 = time.time()
#         print('Integration took: {:.3f}s'.format(t2-t1))
#         if result.y[1][-1]>0:
#             color = 'orange'
#         else:
#             color = 'blue'
#         plt.plot(result.y[0],result.y[1],c=color)
#         stable.append([result.y[0][-1],result.y[1][-1]])
#         
# stable = np.array(stable)
# 
# plt.plot(stable[-1,0],stable[-1,1],'k*',ms=15)
# plt.plot(stable[0,0],stable[0,1],'k*',ms=15)
# plt.xlabel('n')
# plt.ylabel('s')
# =============================================================================
