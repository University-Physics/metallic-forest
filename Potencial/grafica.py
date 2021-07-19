import matplotlib.pyplot as plt
import numpy as np
import re
from fronteras import opening
from fronteras import fractal_dimension
l=input("Escribe...")
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
X=[[] for m in range(4)]
Y=[[] for m in range(4)]
Z=[[] for m in range(4)]
dimension=[[0 for i in range(4)] for i in range(2)]
for m in range(0,4):
    archivo=open(str(m)+l)
    A=archivo.read().split("\n")
    c=opening(str(m)+l)
    #dimension[m]=fractal_dimension(c, max_box_size = None, min_box_size = 1, n_samples = 40, n_offsets = 0, plot = False)
    for line in A:

        C=re.findall("[0-9]*",line)
        if len(C) > 1:
            X[m].append(int(C[0]))
            Y[m].append(int(C[3]))
            Z[m].append(int(C[7]))
    colors = (0,0,0)
    area = np.pi*3

    G=re.findall("[0-9]+",l)
    condicion=G[0]
    T=str((0.001*float(G[0]))**2)
    V=str(float(G[1])*0.1)
    R=str(0.01*float(G[2]))
    fig.suptitle("Growth of dendrites"+r" $T="+T+"\;V="+V+"\;R="+R+"\; N_A=N_B"+"$")
    print(G)
Casos=["0","1","2","3"]
for i in Casos:
    for j in Casos:
        v="condition1T10V1R2I"
        c=opening(i+v+j+"S.txt")
        aux=[0 for i in range(4)]
        aux[int(j)]=fractal_dimension(c, max_box_size = None, min_box_size = 1, n_samples = 40, n_offsets = 0, plot = False)
        dimension[0][int(i)]+=fractal_dimension(c, max_box_size = None, min_box_size = 1, n_samples = 40, n_offsets = 0, plot = False)
        dimension[1][int(i)]=np.std(aux)/2
    dimension[0][int(i)]/=4
ax1.scatter(X[0], Y[0], s=area, c=Z[0], alpha=0.5,label="Contorno 1\n D="+str(round(dimension[0][0],2))+"$\pm$"+str(round(dimension[1][0],2)))
ax1.legend()
ax2.scatter(X[1], Y[1], s=area, c=Z[1], alpha=0.5,label="Contorno 2\n D="+str(round(dimension[0][1],2))+"$\pm$"+str(round(dimension[1][1],2)))
ax2.legend()
ax3.scatter(X[2], Y[2], s=area, c=Z[2], alpha=0.5,label="Contorno 3\n D="+str(round(dimension[0][2],2))+"$\pm$"+str(round(dimension[1][2],2)))
ax3.legend()
ax4.scatter(X[3], Y[3], s=area, c=Z[3], alpha=0.5,label="Contorno 4\n D="+str(round(dimension[0][3],2))+"$\pm$"+str(round(dimension[1][3],2)))
ax4.legend()
for ax in fig.get_axes():
    ax.label_outer()
plt.legend()
plt.savefig("Dentrities.pdf")
plt.show()
    
