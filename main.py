# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import numpy as np

from CLT_Laminate import Laminate
from CLT_Loads import thermal_loads
from CLT_Material import Material
from CLT_Ply import Ply

import matplotlib.pyplot as plt
plt.close('all')

a = Material(E1=181e9,
             E2=10.3e9,
             G12=7.17e9,
             v12=0.28,
             a11=-0.4e-6,
             a22=18e-6,
             name='Gr/Ep',
             verbose=False)

a.addFailureTW(s1t=2.723e3,
               s1c=1.2e3,
               s2t=127,
               s2c=200,
               t12=95.8,
               F12=0)

c = Laminate(stack=[0],
             material=[a],
             thickness=[0.001],
             symmetric=False,
             repeatLeft=1,
             repeatRight=1,
             verbose=True)

F = np.array([250.,  # Nx
              0.,  # Ny
              0.,  # Nxy
              0.,  # Mx
              0.,  # My
              0.  # Mxy
              ])

# F_t = thermal_loads(-155.6, c)
# F = F + F_t
# print(np.array_str(F, precision=2, suppress_small=True))


# x = np.linspace(0, 90)
# y = np.zeros_like(x)
# for i, angle in enumerate(x):
#     print(angle)
#     c = Laminate(stack=[angle],
#                  material=[a],
#                  thickness=[0.005],
#                  symmetric=False,
#                  repeat_left=1,
#                  repeat_right=1,
#                  verbose=True)

#     F = np.array([250.,  # Nx
#                   0.,  # Ny
#                   0.,  # Nxy
#                   0.,  # Mx
#                   0.,  # My
#                   0.  # Mxy
#                   ])
    
#     c.calculateStrain(*F)
#     print(c.lam_strain)
#     y[i] = c.lam_strain[0]

# plt.scatter(x, y)
    
# c.ProgressiveFailureAnalysis(*F, True)



# c.plot_strains(glob=False, idx=0)
# c.plot_curvature()


# If you want to plot just E1 as a function of angles:
# [test] = a.calculateE1Parametric(Angles=linspace(0, 90, 100));
# b.ABD
# figure
# plot(linspace(0, 90, 100), test);


# Example 2.4.2
a = Material(E1=181e9,
             E2=10.3e9,
             G12=7.17e9,
             v12=0.28,
             a11=-0.4e-6,
             a22=18e-6,
             name='Gl/Ep',
             verbose=False)

F = np.array([0.,  # Nx
              0.,  # Ny
              0.,  # Nxy
              1.e6,  # Mx
              0.,  # My
              0.  # Mxy
              ])

x = np.linspace(0, 45)
ex = np.zeros_like(x)
ey = ex.copy()
gxy = ex.copy()
kx = ex.copy()
ky = ex.copy()
kxy = ex.copy()

for i, angle in enumerate(x):
    c = Laminate(stack=[angle,
                        angle + 90,
                        angle + 90,
                        angle,
                        angle,
                        angle + 90,
                        angle + 90,
                        angle,
                        45,
                        -45,
                        -45,
                        45,
                        45,
                        -45,
                        -45,
                        45],
                 material=[a],
                 thickness=[0.001],
                 symmetric=False,
                 repeatLeft=1,
                 repeatRight=1,
                 verbose=True)

    c.calculateStrain(*F)

    ex[i] = c.lam_strain[0]
    ey[i] = c.lam_strain[1]
    gxy[i] = c.lam_strain[2]
    kx[i] = c.lam_strain[3]
    ky[i] = c.lam_strain[4]
    kxy[i] = c.lam_strain[5]


plt.figure()
plt.plot(x, ex)
plt.plot(x, -ey)
plt.plot(x, gxy)

plt.figure()
plt.plot(x, kx)
plt.plot(x, -ky)
plt.plot(x, -kxy)

plt.show()
