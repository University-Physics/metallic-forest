import numpy as np
import matplotlib.pyplot as plt
import re
def read_file(filename):
    archivo=open(filename)
    muestras=[]# 8 81818
    dimension=[[]]
    division=[]
    A=archivo.read().split("\n")
    for line in A:
        C=re.findall("[0-9]+.[0-9]+",line)
        if A.index(line)==0:
            if float(C[0]) not in muestras:
                muestras.append(float(C[0]))
                dimension.append([])
                division.append(0.0)
        if len(line)>=1:
            if float(C[0]) in muestras:
                ind=muestras.index(float(C[0]))
                dimension[ind].append(float(C[1]))
                division[ind]+=1
            else:
                muestras.append(float(C[0]))
                ind=muestras.index(float(C[0]))
                dimension.append([])
                division.append(0.0)
                dimension[ind].append(float(C[1]))
                division[ind]+=1
    definitivos=np.zeros((2,len(muestras)))
    for i in range(len(muestras)):
        definitivos[0,i]=np.mean(dimension[i])
        definitivos[1,i]=np.std(dimension[i])/np.sqrt(division[i])
        print(definitivos[1,i])
    return muestras,definitivos
            
def plot_imbalance(filename, l, rep, outfile, labelphrase):
    x,y =read_file(filename)
    for i in range(len(x)):
        x[i]=1-2/x[i]
    plt.figure()
    plt.style.use("seaborn-white")
    plt.xlim([-1.1,1])
    plt.errorbar(x, y[0,:],yerr=y[1,:],ms=7, color='k', fmt=".b",label=labelphrase)
    plt.plot(x, y[0,:],color='k')
    a0=1.603
    a1=-0.43
    a2=1.78
    X=np.linspace(-1,0.8,num=100)
    Y=a0+a1*X/(a2-X)**2
    plt.plot(X,Y,color="g",label="fit")
    plt.text(-1,1.475,r"$y(x)=a_0+a_1\cdot\frac{x}{(x-a_2)^2}$")
    plt.text(-1,1.44,r"$a_0=1.603\pm 0.004$")
    plt.text(-1,1.42,r"$a_1=-0.43\pm 0.05$")
    plt.text(-1,1.40,r"$a_2=1.78\pm 0.07$")
    plt.xlabel(r'$\frac{N_A-N_B}{N_A+N_B}$')
    plt.ylabel(r'Dimensión fractal')
    plt.legend()
    plt.savefig(outfile+".png")
    for i in range(len(x)):
        with open(outfile+".txt", "a") as text_file:
            text_file.write(str(x[i])+ "\t"+str(y[0,i])+"\n")
    return

plot_imbalance("Results_imbalance100T100V1R.txt", 10, 10, "I4", "kT="+str(0.01)+", R/L=0.01, V=10.0")

