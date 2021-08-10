import matplotlib.pyplot as plt
from leer_archivos import read_file
from leer_archivos import read_Probability
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
    if tipo=="frontera":
        xlabel="Número de frontera"
    plt.figure()
    plt.errorbar(x,y[0,:],yerr=y[1,:],color='black', label=labelphrase)
    plt.scatter(x,y[0,:], color='red')
    plt.xlabel(xlabel)
    plt.ylabel('Dimensión fractal')
    plt.legend()
    plt.savefig(outfile+".png")
    for i in range(len(x)):
        with open(outfile+".txt", "a") as text_file:
            text_file.write(str(x[i])+ "\t"+str(y[0,i])+"\t"+str(y[1,i])+"\n")
    return
def graph_probability(filename):
    I=["0","1","2","3"]
    G=["blue","red","green","orange"]
    labelphrase="contorno "
    plt.figure()
    plt.xlabel("iteración")
    plt.ylabel('Partículas electrodepositadas')
    plt.xlim([0,1080])
    for i in I:
        x,y=read_Probability(i+filename)
        plt.errorbar(x,y[0,:],color=G[int(i)], label=labelphrase+i+" S=2")
        plt.scatter(x,y[0,:], color=G[int(i)])
    plt.legend()
    plt.savefig("Particulas_electrodepositadas2.png")
        
