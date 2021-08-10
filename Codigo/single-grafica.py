import matplotlib.pyplot as plt
import numpy as np
from leer_archivos import read_file
def plot(filename,tipo,outfile, labelphrase):
    x,y =read_file("w"+filename)
    x1,y1=read_file("wo"+filename)
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
            x1[i]=0.02*x1[i]
            xlabel=r'$\frac{\Delta V q_A}{|R_A+R_B|}$'
    plt.figure()
    plt.style.use("seaborn-white")
    plt.errorbar(x,y[0,:],yerr=y[1,:],color='blue', label=labelphrase+" Con I")
    plt.scatter(x,y[0,:], color='black')
    plt.errorbar(x1,y1[0,:],yerr=y1[1,:],color='red', label=labelphrase+" Sin I")
    plt.scatter(x1,y1[0,:], color='black')
    csfont = {'fontname':'Comic Sans MS'}
    
    plt.text(0.01,1.25,"Régimen\n dispersivo",**csfont,fontsize=14)
    plt.xlim([0,0.4])
    plt.xlabel(xlabel)
    plt.ylabel('Dimensión fractal')
    plt.legend()
    plt.savefig(outfile+".png")
plot("Results_potential1T1R2I.txt","potential", "V2", "kT="+str(0.000001)+", R/L=0.001, I=0")

#To plot the I's files we define the expression:
#xe=np.linspace(-1,0.8,100)
#a=1.95
#b=-0.58
#c=1.63
#ye=[a+b/(c-x) for x in xe]
 
#I1
# plt.text(-0.76,1.57,r"$a=(182\pm 1)\cdot 10^{-2}$")
# plt.text(-0.76,1.54,r"$b=-(16 \pm 1)\cdot 10^{-2}$")
# plt.text(-0.76,1.51,r"$c=(106 \pm 5)\cdot 10^{-2}$")

#I2
#plt.text(-0.76,1.6,r"$y(x)=a+b/(c-x)$")
#plt.text(-0.76,1.57,r"$a=(182\pm 1)\cdot 10^{-2}$")
#plt.text(-0.76,1.54,r"$b=-(16 \pm 1)\cdot 10^{-2}$")
#plt.text(-0.76,1.51,r"$c=(105 \pm 5)\cdot 10^{-2}$")

#I3
#plt.text(-0.76,1.6,r"$y(x)=a+b/(c-x)$")
#plt.text(-0.76,1.57,r"$a=(185\pm 3)\cdot 10^{-2}$")
#plt.text(-0.76,1.54,r"$b=-(35 \pm 6)\cdot 10^{-2}$")
#plt.text(-0.76,1.51,r"$c=(140 \pm 9)\cdot 10^{-2}$")

#I4
#plt.text(-0.76,1.6,r"$y(x)=a+b/(c-x)$")
#plt.text(-0.76,1.57,r"$a=(195\pm 5)\cdot 10^{-2}$")
#plt.text(-0.76,1.54,r"$b=-(58 \pm 1)\cdot 10^{-1}$")
#plt.text(-0.76,1.51,r"$c=(16 \pm 1)\cdot 10^{-1}$")
