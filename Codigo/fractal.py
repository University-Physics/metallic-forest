import matplotlib.pyplot as plt
import numpy as np
import re
from fronteras import opening
from Computar_fractal import fractal_dimension
l=input("Escribe...")
X=[]
Y=[]
Z=[]
dimension=0
archivo=open(l)
A=archivo.read().split("\n")
c=opening(l)
D=re.findall("[0-9]+",l)
if int(D[0])==0:
    dimension=1.656
    error=0.001
if int(D[0])==1:
    dimension=1.597
    error=0.001
if int(D[0])==2:
    dimension=1.629
    error=0.002
if int(D[0])==3:
    dimension=1.596
    error=0.001
#dimension=fractal_dimension(c, max_box_size = None, min_box_size = 1, n_samples = 40, n_offsets = 0, plot = False)
for line in A:

    C=re.findall("[0-9]+",line)
    if len(C) > 1:
        print(C)
        if int(C[2])==1:
            X.append(int(C[0]))
            Y.append(int(C[1]))
            Z.append(int(C[2]))
    G=re.findall("[0-9]+",l)
    T=str(0.001*float(G[1]))
    V=str(float(G[2])*0.1)
    R=str(0.01*float(G[3]))
    I=str(float(G[4]))
    plt.title("Growth of dendrites"+r" $T="+T+"\;V="+V+"\;R="+R+"\;I="+I+"$",fontsize=15)
    plt.style.use("classic")
colors = (0,1)

area = np.pi*3
plt.xlim([0,200])
plt.ylim([0.,200])
plt.scatter(X, Y, s=area, c="brown", alpha=0.5,label=r"$D="+str(dimension)+"\pm"+str(error)+"$")
plt.legend(fontsize=20)
plt.savefig("Dendrite.png")
plt.show()
    
