import numpy as np
import glob
from scipy import optimize
from utils import *
from tqdm import tqdm

folder = "./SimplexAlgorithm/data/netlib/*.mps"

good_filenames = glob.glob(folder)

primal_files = []
revised_files = []

for file in tqdm(good_filenames):
    try:
        A, b, c = load_mps(file)
    except:
        continue
    m, n = A.shape
    if (m < 1000 and n < 1000):
        try:
            res = optimize.linprog(c, A_ub = A, b_ub = b, method='revised simplex')
            if res.success:
                revised_files.append(file)
        except:
            pass

        try:
            res = optimize.linprog(c, A_ub = A, b_ub = b, method='simplex')
            if res.success:
                primal_files.append(file)
        except:
            pass

print("Primal:", primal_files)
print("Revised:",revised_files)