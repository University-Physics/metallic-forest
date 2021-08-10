import matplotlib.pyplot as plt
from leer_archivos import read_file
def plot(filename,tipo,outfile, labelphrase):
    x,y =read_file(filename)
    if tipo=="imbalance":
        for i in range(len(x)):
            x[i]=1-2/x[i]
            xlabel=r'$\frac{N_A-N_B}{N_A+N_B}$'
    if tipo=="temperature":
        for i in range(len(x)):
            x[i]=(0.001*x[i])**2
            xlabel=r'$kT$'
    if tipo=="radii":
        for i in range(len(x)):
            x[i]=0.001*x[i]
            xlabel=r'$r/L$'
    if tipo=="potential":
        for i in range(len(x)):
            x[i]=0.02*x[i]
            xlabel=r'$\frac{\Delta V q_A}{|R_A+R_B|}$'
    plt.figure()
    plt.errorbar(x,y[0,:],yerr=y[1,:],color='black', label=labelphrase)
    plt.scatter(x,y[0,:], color='red')
    plt.xlabel(xlabel)
    plt.ylabel('Dimensión fractal')
    plt.legend()
    plt.savefig(outfile+".png")
    for i in range(len(x)):
        with open(outfile+".txt", "a") as text_file:
            text_file.write(str(x[i])+ "\t"+str(y[0,i])+"\n")
    return
