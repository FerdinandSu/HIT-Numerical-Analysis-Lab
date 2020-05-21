import numpy as np
from numpy.linalg import solve

x=[19,25,31,38,44]
y=[19.0,32.3,49.0,73.3,97.8]

phi1=[i*i for i in x]
phi0=[1]*5

phi00=5
phi11=sum([i**2 for i in phi1])
phi01=sum(phi1)
yphi0=sum(y)
yphi1=sum([y[i]*phi1[i] for i in range(0,5)])

print(phi00)
print(phi01)
print(phi11)
print(yphi0)
print(yphi1)

a=np.mat([[phi00,phi01],[phi01,phi11]])#系数矩阵
b=np.mat([yphi0,yphi1]).T    #常数项列矩阵
x,y=solve(a,b)        #方程组的解

print(x,"\n",y)