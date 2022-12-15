import numpy as np
from scipy import optimize
import primal_simplex
import revised_simplex
from utils import *
from tqdm import tqdm
from timeit import default_timer as timer

# TESTS: Afiro, Adlittle, agg2, capri

primal_test_files = ['./SimplexAlgorithm/data/netlib\\adlittle.mps', './SimplexAlgorithm/data/netlib\\afiro.mps', './SimplexAlgorithm/data/netlib\\agg2.mps', './SimplexAlgorithm/data/netlib\\beaconfd.mps', './SimplexAlgorithm/data/netlib\\capri.mps', './SimplexAlgorithm/data/netlib\\lotfi.mps']
revised_test_files = ['./SimplexAlgorithm/data/netlib\\adlittle.mps', './SimplexAlgorithm/data/netlib\\afiro.mps', './SimplexAlgorithm/data/netlib\\agg.mps', './SimplexAlgorithm/data/netlib\\agg2.mps', './SimplexAlgorithm/data/netlib\\beaconfd.mps', './SimplexAlgorithm/data/netlib\\capri.mps']

ptimes = open("./SimplexAlgorithm/Data/ptimes.txt", "w")
rtimes = open("./SimplexAlgorithm/Data/rtimes.txt", "w")
ptimes.write("Scipy time,My time\n")
rtimes.write("Scipy time,My time\n")
numiter = 5

# Test Primal Algorithm
for file in tqdm(primal_test_files):
    A, b, c = load_mps(file)
    avg_sci_time = 0
    res = None
    for i in range(numiter):
        start = timer()
        res = optimize.linprog(c, A_ub = A, b_ub = b, method='simplex')
        end = timer()
        avg_sci_time += end - start
    avg_sci_time /= numiter
    sci_val = res.fun
    avg_my_time = 0
    
    myval = 0
    for i in range(numiter):
        start = timer()
        _, myval = primal_simplex.primal_simplex(A, b, c)
        end = timer()
        avg_my_time += end - start
    avg_my_time /= numiter
    ptimes.write(np.format_float_positional(avg_sci_time, precision=3) + "," + np.format_float_positional(avg_my_time, precision=3) + ",\n")

for file in tqdm(revised_test_files):
    A, b, c = load_mps(file)
    avg_sci_time = 0
    res = None
    for i in range(numiter):
        start = timer()
        res = optimize.linprog(c, A_ub = A, b_ub = b, method='revised simplex')
        end = timer()
        avg_sci_time += end - start
    avg_sci_time /= numiter
    sci_val = res.fun
    avg_my_time = 0
    
    myval = 0
    for i in range(numiter):
        start = timer()
        _, myval = revised_simplex.revised_simplex(A, b, c)
        end = timer()
        avg_my_time += end - start
    avg_my_time /= numiter
    rtimes.write(np.format_float_positional(avg_sci_time, precision=3) + "," + np.format_float_positional(avg_my_time, precision=3) + ",\n")
ptimes.close()
rtimes.close()

print("FINISHED PROFILING")
# # 
# A, b, c = load_mps("./SimplexAlgorithm/data/netlib/afiro.mps")
# _, sol = revised_simplex.revised_simplex(A, b, c)
# print(sol)

# res = optimize.linprog(c, A_ub = A, b_ub = b, method='revised simplex')
# print(res.fun, res.success)

