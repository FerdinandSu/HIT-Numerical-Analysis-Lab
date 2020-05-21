import numpy as np
from math import *
from numpy.linalg import solve

phi0 = 1
phi01 = 1/2
phi02 = 1/3
phi1 = 1/3
phi12 = 1/4
phi2 = 1/5

yphi0 = 2/pi
yphi1 = 1/pi
yphi2 = 1/pi-4/pi**3


a = np.mat([
    [phi0, phi01, phi02],
    [phi01, phi1, phi12],
    [phi02,phi12,phi2]])  # 系数矩阵
b = np.mat([yphi0, yphi1,yphi2]).T  # 常数项列矩阵
x = solve(a, b)  # 方程组的解

print(x)
