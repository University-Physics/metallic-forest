import numpy as np
import matplotlib.pyplot as plt
import re
from funciones import fractal_dimension
from funciones import opening

name=input("Hello, say the name of the file: ")
dimension=fractal_dimension(opening(name), max_box_size = None, min_box_size = 1, n_samples = 40, n_offsets = 0, plot = True)
