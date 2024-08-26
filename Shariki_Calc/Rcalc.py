from genericpath import samefile
from inspect import signature
import numpy as np
import numpy.random as rand
from numpy import sin,pi
import matplotlib.pyplot as plt
import subprocess
import math

def CalcR(Sigma, Density, G, a):
   R = pow((2*Sigma/(Density*G)), 0.5)/a
   return R

def CalcDeltaR(Sigma, Density, G, a):
   R = pow((2*Sigma/(Density*G)), 0.5)/a
   return R

Sigma = 28.63/1000
Density = 861.12
G = 9.81

a = 0.59
R = 1000 * CalcR(Sigma, Density, G, a)
print('a = ', a, 'R =', R, 'mm')
print('Если размер форсунки больше, чем этот радиус, то будет наблюдаться неустойчивость Тейлора. В таком случае размеры шариков будут случайны')
print()

a = 1
R = 1000 * CalcR(Sigma, Density, G, a)
print('a = ', a, 'R max =', R, 'mm')
a = 1.265
R = 1000 * CalcR(Sigma, Density, G, a)
print('a = ', a, 'R min =', R, 'mm')
print('Если размер форсунки меньше максимального, но больше минимального радиуса, то можно наблюдать примерно устойчивый отрыв')
print()

a = 1.265
R = 1000 * CalcR(Sigma, Density, G, a)
print('a = ', a, 'R =', R, 'mm')
print('Если размер форсунки меньше, чем этот радиус, то отрыв будет возникать только в процессе развития неустойчивость, а значит радиус капли будет варироваться')
print()

Radius = 4 / 1000
Kek = Density * G * 3.14 * 4 * pow(Radius, 3) / 3
Lol = 2 * 3.14 * Sigma * Radius/2
print(Kek)
print(Lol)
# dSigma = 
# dDensity = 
# dG = 
# da = 
# DeltaR = CalcDeltaR(Sigma, Density, G, a, dSigma, dDensity, dG, da)
# print('dR =', DeltaR, 'mm')