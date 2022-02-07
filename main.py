# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import numpy as np

from CLT_Laminate import Laminate
from CLT_Loads import thermal_loads
from CLT_Material import Material
from CLT_Ply import Ply

a = Material(e1=25e6,
             e2=1.5e6,
             g12=1.9e6,
             v12=0.3,
             a11=-0.5e-6,
             a22=15e-6,
             name='Gr/Ep',
             verbose=False)

b = Material(e1=8e6,
             e2=2.3e6,
             g12=1.1e6,
             v12=0.28,
             a11=3.7e-6,
             a22=14e-6,
             name='Gl/Ep',
             verbose=False)

a.addFailureTW(s1t=2.723e3,
               s1c=1.2e3,
               s2t=127,
               s2c=200,
               t12=95.8,
               F12=0)

c = Laminate(stack=[0, 30, 90, -30],
             material=[a, b, a, b],
             thickness=[0.005],
             symmetric=False,
             repeat_left=1,
             repeat_right=1,
             verbose=True)

F = np.array([520,
              377,
              64.4,
              -4.0,
              0.22,
              -0.0854])
F_t = thermal_loads(-275, c)
F = F + F_t
# print(F)
# print(c.ABD)
print(c.abd)
print(c.calculateStrain(*F))
c.plyStrains()
# b.calculateBuckling(Lx=100, ...
# Ly = 100, ...
# k = 0.5, ...
# BC = 'S4', ...
# DisplayPlot = true, ...
# DisplayTable = true);
# print(
#     b.calculateStrain(Nx=10,
#                       Ny=0,
#                       Nxy=0,
#                       Mx=0,
#                       My=0,
#                       Mxy=0)
# )
# b.ProgressiveFailureAnalysis(Nx=32,
#                              Ny=0,
#                              Nxy=0,
#                              Mx=1000,
#                              My=0,
#                              Mxy=0,
#                              verbose=True)

# [a, ~, c, ~] = b.calculateMacroBehavior();

# If you want to plot just E1 as a function of angles:
# [test] = a.calculateE1Parametric(Angles=linspace(0, 90, 100));
# b.ABD
# figure
# plot(linspace(0, 90, 100), test);


# def print_hi(name):
#     # Use a breakpoint in the code line below to debug your script.
#     print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.
#
#
# # Press the green button in the gutter to run the script.
# if __name__ == '__main__':
#     print_hi('PyCharm')

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
