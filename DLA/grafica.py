import matplotlib.pyplot as plt
import numpy as np
import re
archivo=open("data.txt")
A=archivo.read().split("\n")
X=[]
Y=[]
for line in A:
    if len(re.findall("\s0",line))<1:
        C=re.findall("[0-9]*",line)
        if len(C) > 1:
            X.append(int(C[0]))
            Y.append(int(C[2]))
colors = (0,0,0)
area = np.pi*3
plt.scatter(X, Y, s=area, c=colors, alpha=0.5)
plt.title("DLA")
plt.xlabel('x')
plt.ylabel('y')
plt.ylim([0,100])
plt.show()
plt.savefig("DLA.png")
