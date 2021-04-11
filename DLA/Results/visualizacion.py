# 1 "./grafica.py"
# 1 "<built-in>"
# 1 "<command-line>"
# 31 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 32 "<command-line>" 2
# 1 "./grafica.py"
import matplotlib.pyplot as plt
import numpy as np
import re
for i, n in enumerate(["100.000","250.000","500.000"]):
    archivo=open(n+".txt")
    A=archivo.read().split("\n")
    X=[]
    Y=[]
    Z=[]
    for line in A:
        line=line.strip()
        C=re.findall("[0-9]*",line)
        
        if len(C) > 1:
            X.append(int(C[0]))
            Y.append(int(C[2]))
            Z.append(int(C[4]))
    colors = (0,0,0)
    area = np.pi*3
    plt.set_cmap("Greys")
    plt.scatter(X, Y, s=area, c=Z, alpha=0.5)
    plt.title("DLA")
    plt.xlabel('x')
    plt.ylabel('y')
    plt.ylim([200,400])
    plt.xlim([320,520])
    plt.show()
plt.savefig("DLA.png")
