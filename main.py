# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import numpy as np

from CLT_Laminate import Laminate
from CLT_Loads import thermal_loads
from CLT_Material import Material
from CLT_Ply import Ply

a = Material(e1=20.01e6,
             e2=1.301e6,
             g12=1.001e6,
             v12=0.3,
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

c = Laminate(stack=[0, 45],
             material=[a],
             thickness=[0.005],
             symmetric=False,
             repeat_left=1,
             repeat_right=1,
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
print(np.array_str(F, precision=2, suppress_small=True))

c.calculateStrain(*F)
# c.ProgressiveFailureAnalysis(*F, True)



c.plot_strains(glob=False, idx=0)
c.plot_curvature()


# If you want to plot just E1 as a function of angles:
# [test] = a.calculateE1Parametric(Angles=linspace(0, 90, 100));
# b.ABD
# figure
# plot(linspace(0, 90, 100), test);


