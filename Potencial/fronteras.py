import numpy as np
import matplotlib.pyplot as plt
import re
from funciones import fractal_dimension
from funciones import opening
fron funciones import Readf

Casos=["0","1","2","3"]
equis=[]
for i in Casos:
    dimension=0
    for j in Casos:
        v="condition1T10V1R2I"
        c=opening(i+v+j+"S.txt")
        dimension+=fractal_dimension(c, max_box_size = None, min_box_size = 1, n_samples = 40, n_offsets = 0, plot = False)
    dimension/=4
    equis.append(dimension)
    print(i,dimension)


