import numpy as np
import matplotlib.pyplot as plt
import re

def opening(textname):
    archivo=open(textname)
    A=archivo.read().split("\n")
    X=[]
    Y=[]
    Z=[]
    data=np.zeros(shape=(201, 201))
    for line in A:

        C=re.findall("[0-9]+",line)
        if len(C) > 1:
            X.append(int(C[0]))
            Y.append(int(C[1]))
            Z.append(int(C[2]))
            data[X[-1], Y[-1]]=Z[-1]
    return data

def fractal_dimension(array, max_box_size = None, min_box_size = 1, n_samples = 40, n_offsets = 0, plot = False):
    """Calculates the fractal dimension of a 3D numpy array.
    
    Args:
        array (np.ndarray): The array to calculate the fractal dimension of.
        max_box_size (int): The largest box size, given as the power of 2 so that
                            2**max_box_size gives the sidelength of the largest box.                     
        min_box_size (int): The smallest box size, given as the power of 2 so that
                            2**min_box_size gives the sidelength of the smallest box.
                            Default value 1.
        n_samples (int): number of scales to measure over.
        n_offsets (int): number of offsets to search over to find the smallest set N(s) to
                       cover  all voxels>0.
        plot (bool): set to true to see the analytical plot of a calculation.
                            
        
    """
    #determine the scales to measure on
    if max_box_size == None:
        #default max size is the largest power of 2 that fits in the smallest dimension of the array:
        max_box_size = int(np.floor(np.log2(np.min(array.shape))))
    scales = np.floor(np.logspace(max_box_size,min_box_size, num = n_samples, base =2 ))
    scales = np.unique(scales) #remove duplicates that could occur as a result of the floor
    
    #get the locations of all non-zero pixels
    locs = np.where(array > 0)
    voxels = np.array([(x,y) for x,y in zip(*locs)])
    
    #count the minimum amount of boxes touched
    Ns = []
    #loop over all scales
    for scale in scales:
        touched = []
        if n_offsets == 0:
            offsets = [0]
        else:
            offsets = np.linspace(0, scale, n_offsets)
        #search over all offsets
        for offset in offsets:
            bin_edges = [np.arange(0, i, scale) for i in array.shape]
            bin_edges = [np.hstack([0-offset,x + offset]) for x in bin_edges]
            H1, e = np.histogramdd(voxels, bins = bin_edges)
            touched.append(np.sum(H1>0))
        Ns.append(touched)
    Ns = np.array(Ns)
    
    #From all sets N found, keep the smallest one at each scale
    Ns = Ns.min(axis=1)
   
    
    
    #Only keep scales at which Ns changed
    scales  = np.array([np.min(scales[Ns == x]) for x in np.unique(Ns)])
    
    
    Ns = np.unique(Ns)
    Ns = Ns[Ns > 0]
    scales = scales[:len(Ns)]
    #perform fit
    coeffs = np.polyfit(np.log(1/scales), np.log(Ns),1)
    
    #make plot
    if plot:
        fig, ax = plt.subplots(figsize = (8,6))
        ax.scatter(np.log(1/scales), np.log(np.unique(Ns)), c = "teal", label = "Measured ratios")
        ax.set_ylabel("$\log N(\epsilon)$")
        ax.set_xlabel("$\log 1/ \epsilon$")
        fitted_y_vals = np.polyval(coeffs, np.log(1/scales))
        ax.plot(np.log(1/scales), fitted_y_vals, "k--", label = f"Fit: {np.round(coeffs[0],3)}X+{coeffs[1]}")
        ax.legend()
        ax.show();
    print(coeffs[0])
    return(coeffs[0])

def read_file(filename):
    archivo=open(filename)
    muestras=[]
    dimension=[]
    division=[]
    A=archivo.read().split("\n")
    for line in A:
        C=re.findall("[0-9]+.[0-9]+",line)
        if A.index(line)==0:
            if float(C[0]) not in muestras:
                muestras.append(float(C[0]))
                dimension.append(0.0)
                division.append(0.0)
        if len(line)>=1:
            if float(C[0]) in muestras:
                ind=muestras.index(float(C[0]))
                dimension[ind]+=float(C[1])
                division[ind]+=1
            else:
                muestras.append(float(C[0]))
                ind=muestras.index(float(C[0]))
                dimension.append(0.0)
                division.append(0.0)
                dimension[ind]+=float(C[1])
                division[ind]+=1
    for i in range(len(muestras)):
        dimension[i]/=division[i]
        print(muestras[i],dimension[i])
    return muestras,dimension
            

def read_imbalance(T, V, R):
    data=np.zeros((1,2))
    outfile="Results_imbalance"+str(int(T))+"T"+str(int(V))+"V"+str(int(R))+"R.txt"
    for j in range(10):
        for i in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]:
            filename="data/Out"+str(int(T))+"T"+str(int(V))+"V"+str(int(R))+"R"+str(int(i))+"I"+str(int(j+1))+"S.txt"
            data[0,0]=i
            data[0,1]=fractal_dimension(opening(filename), max_box_size = None, min_box_size = 1, n_samples = 100, n_offsets = 0, plot = False)
            with open(outfile, "a") as text_file:
                text_file.write(str(data[0,0])+ "\t"+str(data[0,1])+"\n")
    return

def read_temperature(V, R, I):
    data=np.zeros((1,2))
    outfile="Results_temperature"+str(int(V))+"V"+str(int(R))+"R"+str(int(I))+"I.txt"
    for j in range(10):
        for i in [1, 2, 5, 7, 10, 15, 20, 30, 40, 50, 70, 100]:
            filename="data/Out"+str(int(i))+"T"+str(int(V))+"V"+str(int(R))+"R"+str(int(I))+"I"+str(int(j+1))+"S.txt"
            data[0,0]=i
            data[0,1]=fractal_dimension(opening(filename), max_box_size = None, min_box_size = 1, n_samples = 100, n_offsets = 0, plot = False)
            with open(outfile, "a") as text_file:
                text_file.write(str(data[0,0])+ "\t"+str(data[0,1])+"\n")
    return

def read_radii(T, V, I):
    data=np.zeros((1,2))
    outfile="Results_radii"+str(int(T))+"T"+str(int(V))+"V"+str(int(I))+"I.txt"
    for j in range(10):
        for i in [1, 2, 3, 4, 5]:
            filename="data/Out"+str(int(T))+"T"+str(int(V))+"V"+str(int(i))+"R"+str(int(I))+"I"+str(int(j+1))+"S.txt"
            data[0,0]=i
            data[0,1]=fractal_dimension(opening(filename), max_box_size = None, min_box_size = 1, n_samples = 100, n_offsets = 0, plot = False)
            with open(outfile, "a") as text_file:
                text_file.write(str(data[0,0])+ "\t"+str(data[0,1])+"\n")
    return

def read_potential(T, R, I):
    data=np.zeros((1,2))
    outfile="Results_potential"+str(int(T))+"T"+str(int(R))+"R"+str(int(I))+"I.txt"
    for j in range(10):
        for i in [1, 2, 5, 7, 10, 15, 20, 30, 40, 50, 70, 100, 200, 250, 300, 350, 400, 500, 600, 700]:
            filename="data/Out"+str(int(T))+"T"+str(int(i))+"V"+str(int(R))+"R"+str(int(I))+"I"+str(int(j+1))+"S.txt"
            data[0,0]=i
            data[0,1]=fractal_dimension(opening(filename), max_box_size = None, min_box_size = 1, n_samples = 100, n_offsets = 0, plot = False)
            with open(outfile, "a") as text_file:
                text_file.write(str(data[0,0])+ "\t"+str(data[0,1])+"\n")
    return

def plot_imbalance(filename, l, rep, outfile, labelphrase):
    x,y =read_file(filename)
    plt.figure()
    plt.scatter(x,y, color='k', label=labelphrase)
    plt.plot(x,y, color='k')
    plt.xlabel(r'$\frac{N_A-N_B}{N_A+N_B}$')
    plt.ylabel(r'Dimensión fractal')
    plt.legend()
    plt.savefig(outfile+".png")
    for i in range(len(x)):
        with open(outfile+".txt", "a") as text_file:
            text_file.write(str(x[i])+ "\t"+str(y[i])+"\n")
    return

def plot_temperature(filename, l, rep, outfile, labelphrase):
    x,y =read_file(filename)
    plt.figure()
    plt.scatter(x, y, color='k', label=labelphrase)
    plt.plot(x, y, color='k')
    plt.xlabel(r'$kT$')
    plt.ylabel(r'Dimensión fractal')
    plt.legend()
    plt.savefig(outfile+".png")
    for i in range(len(x)):
        with open(outfile+".txt", "a") as text_file:
            text_file.write(str(x[i])+ "\t"+str(y[i])+"\n")
        
    return

def plot_radii(filename, l, rep, outfile, labelphrase):
    x,y =read_file(filename)
    plt.figure()
    plt.scatter(x, y, color='k', label=labelphrase)
    plt.plot(x, y, color='k')
    plt.xlabel(r'$r/L$')
    plt.ylabel(r'Dimensión fractal')
    plt.legend()
    plt.savefig(outfile+".png")
    for i in range(len(x)):
        with open(outfile+".txt", "a") as text_file:
            text_file.write(str(x[i])+ "\t"+str(y[i])+"\n")
    return

def plot_potential(filename, l, rep, outfile, labelphrase):
    x,y =read_file(filename)
    plt.figure()
    plt.scatter(x, y, color='k', label=labelphrase)
    plt.plot(x, y, color='k')
    plt.xlabel(r'$\frac{\Delta V q_A}{|R_A+R_B|}$')
    plt.ylabel(r'Dimensión fractal')
    plt.legend()
    plt.savefig(outfile+".png")
    for i in range(len(x)):
        with open(outfile+".txt", "a") as text_file:
            text_file.write(str(x[i])+ "\t"+str(y[i])+"\n")
    return

read_imbalance(1, 1, 1)
plot_imbalance("Results_imbalance1T1V1R.txt", 10, 10, "I1", "kT=0.001, R/L=0.01, V=0.1")
read_imbalance(100, 1, 1)
plot_imbalance("Results_imbalance100T1V1R.txt", 10, 10, "I2", "kT=0.1, R/L=0.01, V=0.1")
read_imbalance(1, 100, 1)
plot_imbalance("Results_imbalance1T100V1R.txt", 10, 10, "I3", "kT=0.001, R/L=0.01, V=10.0")
read_imbalance(100, 100, 1)
plot_imbalance("Results_imbalance100T100V1R.txt", 10, 10, "I4", "kT=0.1, R/L=0.01, V=10.0")
read_radii(1, 1, 2)
plot_radii("Results_radii1T1V2I.txt", 5, 10, "R1", "kT=0.001, V=0.1, I=0") #5  #10 re
read_radii(100, 1, 2)
plot_radii("Results_radii100T1V2I.txt", 5, 10, "R2", "kT=0.1, V=0.1, I=0")
read_radii(1, 100, 2)
plot_radii("Results_radii1T100V2I.txt", 5, 10, "R3", "kT=0.001, V=10.0, I=0")
read_radii(100, 100, 2)
plot_radii("Results_radii100T100V2I.txt", 5, 10, "R4", "kT=0.1, V=10.0, I=0")
read_temperature(1, 1, 2)
plot_temperature("Results_temperature1V1R2I.txt", 12, 10, "T1", "V=0.1, R/L=0.01, I=0")
read_temperature(100, 1, 2)
plot_temperature("Results_temperature100V1R2I.txt", 12, 10, "T2", "V=10.0, R/L=0.01, I=0")
read_potential(1, 1, 2)
plot_potential("Results_temperature1V1R2I.txt", 20, 10, "V1", "kT=0.001, R/L=0.01, I=0")
read_potential(100, 1, 2)
plot_potential("Results_temperature1V1R2I.txt", 20, 10, "V2", "kT=0.1, R/L=0.01, I=0")


