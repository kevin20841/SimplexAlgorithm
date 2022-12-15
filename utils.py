from pysmps import smps_loader as smps
import numpy as np
# Implement Utilities:


# TODO: Detect and find dual
# 

def load_mps(path):
    name, objective_name, row_names, col_names, col_types, types, c, A, rhs_names, rhs, bnd_names, bnd= smps.load_mps(path)
    b = np.array([], dtype = np.float64)
    for key in rhs_names:
        b = np.concatenate((b, np.array(rhs[key])))

    for i in range(len(types)):
        if types[i] =='E':
            A = np.vstack((A, -A[i]))
            b = np.append(b, -b[i])
 
        elif types[i] == "G":
            A[i] *= -1
            b[i] *= -1

    return A, b, c
# TODO: Constraint generation, converting between dual, MPS file handling, using pysmps
# https://people.sc.fsu.edu/~jburkardt/datasets/mps/adlittle.mps
# res = smps.load_mps("./data/adlittle.mps")
