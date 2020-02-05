from pylab import *

import os
import subprocess

wd = os.getcwd()
os.chdir(wd)
subprocess.call(["g++", "mes_impl.cpp"])
subprocess.call(["./a.out", "50"])

X = []
Y = []

f = open("./result.txt", "r")

for line in f:
    data = line.split()
    X.append(float(data[0]))
    Y.append(float(data[1]))

plt.plot(X, Y)

plt.xlabel('x')

plt.ylabel('y')

plt.title('Function')

plt.show()