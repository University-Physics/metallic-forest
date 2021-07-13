import numpy as np
import matplotlib.pyplot as plt
import re
def Readf(name,rows):
    file=open(name,"r")
    cols=[[] for i in range(0,rows)]
    data=np.zeros(shape=(201, 201))
    for line in file:
        aux=line.strip()
        F=re.findall(r"[0-9]+",aux)
        if len(F) > 1:
            for i in range(0,rows):
                try:
                    cols[i].append(float(F[i]))
                    data[cols[0][-1],cols[1][-1]]=cols[2][-1]
                except:
                    continue
    return data
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
    return coeffs[0]
Casos=["0","1","2","3"]
equis=[]
for i in Casos:
    dimension=0
    for j in Casos:
        v="condition1T10V1R2I"
        c=opening(i+v+j+"S.txt")
        dimension+=fractal_dimension(c, max_box_size = None, min_box_size = 1, n_samples = 40, n_offsets = 0, plot = False)
    dimension/=4
    print(i,dimension)


