import numpy as np
from primal_simplex import *


A = np.array([[3, 1], [1, 2], [-2, 2]])

b = np.array([180, 100, 40])
c = np.array([4, 12])

_, sol = primal_simplex(A, b, c)
print(sol)


