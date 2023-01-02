import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

plt.close('all')

n = np.zeros(8)

a = 0.1
b = 0
d = 5
g = 1
gs = 0.1
u = 0.7

#  +-----------+-----------+
#  |           |           |
#  |           |           |
#  |     1     |     2     |
#  |     g+    |           |
#  |           |           |
#  +-----------+-----------+
#  |           |           |
#  |           |           |           |
#  |     3     |     4     |
#  |     g-    |           |
#  |           |           |
#  +-----------+-----------+

def ndot(t, n):
    np1, np2, np3, np4 = n[:4]
    nn1, nn2, nn3, nn4 = n[4:]
    
    np1d = a*(np2 + np3 - 4*np1) + b*nn1 - d*np1*nn1 + 0.5*g*(1 + u*(np1 - nn1)/(np1 + nn1)) + gs
    np2d = a*(np1 + np4 - 4*np2) + b*nn2 - d*np2*nn2 + 0.5*g*(1 + u*(np2 - nn2)/(np2 + nn2))
    np3d = a*(np1 + np4 - 4*np3) + b*nn3 - d*np3*nn3 + 0.5*g*(1 + u*(np3 - nn3)/(np3 + nn3))
    np4d = a*(np2 + np3 - 4*np4) + b*nn4 - d*np4*nn4 + 0.5*g*(1 + u*(np4 - nn4)/(np4 + nn4))
    
    nn1d = a*(nn2 + nn3 - 4*nn1) + b*np1 - d*np1*nn1 + 0.5*g*(1 - u*(np1 - nn1)/(np1 + nn1))
    nn2d = a*(nn1 + nn4 - 4*nn2) + b*np2 - d*np2*nn2 + 0.5*g*(1 - u*(np2 - nn2)/(np2 + nn2))
    nn3d = a*(nn1 + nn4 - 4*nn3) + b*np3 - d*np3*nn3 + 0.5*g*(1 - u*(np3 - nn3)/(np3 + nn3)) + gs
    nn4d = a*(nn2 + nn3 - 4*nn4) + b*np4 - d*np4*nn4 + 0.5*g*(1 - u*(np4 - nn4)/(np4 + nn4))
    
    return np.array([np1d, np2d, np3d, np4d, nn1d, nn2d, nn3d, nn4d])

# def ndot(t, n):
#     np1, np2, np3, np4 = n[:4]
#     nn1, nn2, nn3, nn4 = n[4:]
#     c1 = np1 - nn1
#     c2 = np2 - nn2
#     c3 = np3 - nn3
#     c4 = np4 - nn4
    
#     np1d = a*(np2 + np3 - 4*np1) + b*nn1 - d*np1*nn1 + 0.5*g + gs + u*np1*c1
#     np2d = a*(np1 + np4 - 4*np2) + b*nn2 - d*np2*nn2 + 0.5*g + u*np2*c2
#     np3d = a*(np1 + np4 - 4*np3) + b*nn3 - d*np3*nn3 + 0.5*g + u*np3*c3
#     np4d = a*(np2 + np3 - 4*np4) + b*nn4 - d*np4*nn4 + 0.5*g + u*np4*c4
    
#     nn1d = a*(nn2 + nn3 - 4*nn1) + b*np1 - d*np1*nn1 + 0.5*g - u*nn1*c1
#     nn2d = a*(nn1 + nn4 - 4*nn2) + b*np2 - d*np2*nn2 + 0.5*g - u*nn2*c2
#     nn3d = a*(nn1 + nn4 - 4*nn3) + b*np3 - d*np3*nn3 + 0.5*g + gs - u*nn3*c3
#     nn4d = a*(nn2 + nn3 - 4*nn4) + b*np4 - d*np4*nn4 + 0.5*g- u*nn4*c4
    
#     return np.array([np1d, np2d, np3d, np4d, nn1d, nn2d, nn3d, nn4d])


y0 = np.zeros((8))
# konverguje k "s < 0"
y0[2:4] = 1
y0[4:6] = 1
# konvergujek "s > 0"
# y0[0:2] = 0.01
# y0[6:8] = 0.01
out = solve_ivp(ndot, t_span=(0,1), y0=y0)

plt.close('all')
fig, ax = plt.subplots(2, 1, sharex=True)
for nx in out.y:
    ax[0].plot(out.t, nx)
ax[1].semilogy(out.t, np.sum(out.y, axis=0))
    
n_pos = out.y[:4,-1].reshape((2,2))
n_neg = out.y[4:,-1].reshape((2,2))
circ = (n_pos - n_neg)
print(circ)
scale = abs(circ.flatten()).max()
# scale = 1
fig_22, ax_22 = plt.subplots(1, 1)
im = ax_22.imshow(circ, aspect='auto', cmap='seismic', vmin=-scale, vmax=scale)
fig_22.colorbar(im)
fig_22.tight_layout()