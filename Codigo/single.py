import matplotlib.pyplot as plt
import numpy as np
import re
from fronteras import opening
from fronteras import fractal_dimension
l=input("Escribe...")
X=[]
Y=[]
Z=[]
archivo=open(l)
A=archivo.read().split("\n")
c=opening(l)
    #dimension[m]=fractal_dimension(c, max_box_size = None, min_box_size = 1, n_samples = 40, n_offsets = 0, plot = False)
for line in A:
    C=re.findall("[0-9]*",line)
    if len(C) > 1:
        X.append(int(C[0]))
        Y.append(int(C[3]))
        Z.append(int(C[7]))
plt.figure()
colors = (0,0,0)
area = np.pi*3
plt.scatter(X[0], Y[0], s=area, c=Z[0], alpha=0.5)
plt.legend()
plt.savefig("Dentrities.pdf")
plt.show()
    
