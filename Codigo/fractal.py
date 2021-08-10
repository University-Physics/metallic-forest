import matplotlib.pyplot as plt
import numpy as np
import re
from fronteras import opening
from fronteras import fractal_dimension
l=input("Escribe...")
X=[]
Y=[]
Z=[]
dimension=0
archivo=open(l)
A=archivo.read().split("\n")
c=opening(l)
dimension=fractal_dimension(c, max_box_size = None, min_box_size = 1, n_samples = 40, n_offsets = 0, plot = False)
for line in A:

    C=re.findall("[0-9]*",line)
    if len(C) > 1:
        X.append(int(C[0]))
        Y.append(int(C[3]))
        Z.append(int(C[7]))
    G=re.findall("[0-9]+",l)
    T=str(0.001*float(G[0]))
    V=str(float(G[1])*0.1)
    R=str(0.01*float(G[2]))
    I=str(float(G[3]))
    plt.title("Growth of dendrites"+r" $T="+T+"\;V="+V+"\;R="+R+"\;I="+I+"$")
    print(G)
colors = (0,0,0)
area = np.pi*3
plt.scatter(X, Y, s=area, c=Z, alpha=0.5,label="D="+str(round(dimension,2)))
plt.legend()
plt.savefig("Dendrite.pdf")
plt.show()
    
