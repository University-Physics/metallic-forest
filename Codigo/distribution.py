import numpy as np
import matplotlib.pyplot as plt
import re
n=["0","1","2","3"]
x=[[] for i in range(4)]
y=[[] for i in range(4)]
z=[[] for i in range(4)]
def opening(textname):
    archivo=open(textname)
    A=archivo.read().split("\n")
    X=[]
    Y=[]
    for line in A:

        C=re.findall("[0-9]+",line)
        if len(C) > 1:
            X.append(int(C[0]))
            Y.append(int(C[1]))
    return X,Y
plt.xlabel(r"$iteración$")
plt.ylabel(r"$Particulas\;adheridas$")
for i in n:
    x[int(i)],y[int(i)]=opening(i+"Probability_distribution.txt")
    F=x[int(i)]
    G=y[int(i)]
    plt.plot(F,G,label=i+" Contorno")
plt.legend()
plt.savefig("Electrodeposition_distribution.pdf")
plt.show()
