import numpy as np
import matplotlib.pyplot as plt
import re
from Computar_fractal import fractal_dimension
# Función que lee y convierte los archivos.txt en arreglos
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
# Función auxiliar
def Nombre(T,V,R,I,S):
    return  str(int(T))+"T"+str(int(V))+"V"+str(int(R))+"R"+str(int(I))+"I"+str(int(S))+"S.txt"

# Se hacen los archivos aqui
def generate_txt(Name,A,B,C):
    data=np.zeros((1,3))
    if Name=="imbalance":
        outfile="Results_imbalance"+str(int(A))+"T"+str(int(B))+"V"+str(int(C))+"R.txt"
        Variacion=[1,2,3,4,5,6,7,8,9,10]
    elif Name=="temperature":
        outfile="Results_temperature"+str(int(A))+"V"+str(int(B))+"R"+str(int(C))+"I.txt"
        Variacion=[1,2,5,7,10,15,20,30,40,50,70,100]
    elif Name=="potential":
        Variacion=[1, 2, 3, 4, 5, 6 ,7 , 8 ,9 , 10 ,11, 12 ,13 ,14,15, 16, 17, 18, 19,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100] 
        outfile="Results_potential"+str(int(A))+"T"+str(int(B))+"R"+str(int(C))+"I.txt"
    elif Name=="radii":
        Variacion=[1,2,3,4,5]
        outfile="Results_radii"+str(int(A))+"T"+str(int(B))+"V"+str(int(C))+"I.txt"
    elif Name=="frontera":
        Variacion=[0,1,2,3]
        outfile="Results_frontera"+str(int(A))+"T"+str(int(B))+"V"+str(int(C))+"R2I.txt"
    for j in range(10):
        for i in Variacion:
            if Name=="imbalance":
                filename="data/Out"+Nombre(A,B,C,str(int(i)),str(int(j+1)))
            if Name=="temperature":
                filename="data/Out"+Nombre(str(int(i)),A,B,C,str(int(j+1)))
            if Name=="potential":
                filename="data/Out"+Nombre(A,str(int(i)),B,C,str(int(j+1)))
            if Name=="radii":
                filename="data/Out"+Nombre(A,B,str(int(i)),C,str(int(j+1)))
            if Name=="frontera":
                filename="data/"+str(i)+"Out"+Nombre(A,B,C,"2",str(int(j+1)))
            data[0,0]=i
            data[0,1]=fractal_dimension(opening(filename), max_box_size = None, min_box_size = 1, n_samples = 100, n_offsets = 0, plot = False)[0]
            data[0,2]=fractal_dimension(opening(filename), max_box_size = None, min_box_size = 1, n_samples = 100, n_offsets = 0, plot = False)[1]
            with open(outfile, "a") as text_file:
                text_file.write(str(data[0,0])+ "\t"+str(data[0,1])+"\t"+str(data[0,2])+"\n")
    return

# Se lee el archivo aqui
def read_file(filename):
    archivo=open(filename)
    muestras=[]
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
        
        real=[0 for t in range(100)]
        for k in range(100):
            if(len(dimension[i])>=1):
                for j in range(len(dimension[i])):
                    real[k]+=np.random.normal(dimension[i][j][0],dimension[i][j][1],1)
                real[k]/=len(dimension[i])
        definitivos[0,i]=np.mean(real)
        definitivos[1,i]=np.std(real)/np.sqrt(10)
        print(definitivos)
    return muestras,definitivos
                                      
def read_Probability(filename):
    muestras=[]
    dimension=[]
    division=[]
    for i in range(1,11):
        archivo=open("data/fronteras/P/"+filename+str(i)+"SProbability_distribution.txt")
        A=archivo.read().split("\n")
        for line in A:
            C=re.findall("[0-9]+",line)
            print(C)
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
        definitivos[1,i]=np.std(dimension[i])/np.sqrt(10)
        print(definitivos)
    return muestras,definitivos
                                                  
