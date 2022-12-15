import numpy as np
from scipy.optimize._bglu_dense import LU
from numpy.linalg import solve
def initial_BFS(A, b, tol=1e-8):
    [m, n] = A.shape
    x = np.zeros(n)
    residuals = b - A@x
    A[residuals< 0] = - A[residuals < 0]
    b [ residuals < 0] = -b[residuals<0]
    residuals[residuals<0] *= -1
    nzcon = np.arange(m)
    rowsI, colsI = get_sing_col(A, residuals)
    basis = np.array([], dtype=int)
    mask = np.logical_and(np.isin(rowsI, nzcon), np.logical_not(np.isin(colsI, basis)))
    rowsI = rowsI[mask]    
    colsI = colsI[mask]
    auxrows = nzcon[np.logical_not(np.isin(nzcon, rowsI))]

    auxcols = n + np.arange(len(auxrows))
    basisGuess = np.concatenate((colsI, auxcols))
    basisGuessR = np.concatenate((rowsI, auxrows))
    A = np.hstack((A,  np.zeros((m, len(auxrows)))))
    A[auxrows, auxcols] = 1
    x = np.concatenate((x, np.zeros(len(auxrows))))
    x[basisGuess] = residuals[basisGuessR]/A[basisGuessR, basisGuess]
    c = np.zeros(len(auxrows) + n)
    c[auxcols] = 1

    basis = np.concatenate((basis, basisGuess))
    basis = make_basis_nonsingular(A, basis)

    return A, b, c, basis, x

def make_basis_nonsingular(A, basis):
    [m, n] = A.shape

    a = np.arange(m+n)
    bl = np.zeros(m+n, dtype=bool)
    
    bl[basis] = 1
    op = a[~bl]
    op = op[op < n]
    B = np.zeros((m, n))
    for i in range(len(basis)):
        B[:, i] = A[:, basis[i]]
    #B[:, 0:len(basis)] = A[:, basis]

    # check basis not full rank?
    # keep on iterating until we find a full rank basis.
    r = 0
    for i in range(n):
        new_basis = np.random.permutation(op)[:m-len(basis)]
        for j in range(len(basis)):
            B[:, j] = A[:, basis[j]]
        #B[:, len(basis):] = A[:, new_basis]
        r = np.linalg.matrix_rank(B)
        if r == m:
            break
    return np.concatenate((basis, new_basis))
def get_sing_col(A, b):
    ci = np.nonzero(np.sum(np.abs(A)!=0, axis = 0) == 1)[0]
    col = A[:, ci]
    ri = np.zeros(len(ci), dtype=int)
    nzr, nzc = np.nonzero(col)
    ri[nzc] = nzr

    ssign = A[ri, ci] * b[ri] >=0
    ci = ci[ssign][::-1]
    ri = ri[ssign][::-1]
    uri, fc = np.unique(ri, return_index=True)
    
    return uri, ci[fc]

def phase_1(A, b, tol=1e-8):
    [m, n] = A.shape
    # Construct initial BFS
    A, b, c, basis, x = initial_BFS(A, b, tol=tol)
    
    # Phase 1:
    x, basis = solve_simplex(A, c, x, basis)
    
    if c.dot(x) > tol:
        return "Infeasible!"
    
    # Remove artificial variables from basis
    kr = np.ones(m, dtype=bool)

    for bc in basis[basis >=n]:
        B = A[:, basis]
        bf = np.abs(solve(B, A))
        rows = np.argmax(bf[:, bc])
        ec = np.ones(n, dtype=bool)
        ec[basis[basis < n]] = 0
        eci = np.where(ec)[0]
        
        i = np.argmax(bf[:, :n][rows, ec])
        nbc = eci[i]
        if bf[rows, nbc] <tol:
            kr[rows] = False
        else:
            basis[basis == bc] = nbc
    A = A[kr, :n]
    basis = basis[kr]
    x = x[:n]
    
    return A, x, basis

def revised_simplex(A, b, c, tol=1e-8):
    [m, n] = A.shape
    A2 = np.eye(A.shape[0])

    c = np.concatenate((c, np.zeros((A.shape[0],))))
    A = np.hstack([A, A2])
    # Phase 1
    A, x, basis = phase_1(A, b, tol=tol)
    # Phase 2
    x, basis = solve_simplex(A, c, x, basis)

    return x, c.dot(x)

def solve_simplex(A, c, x, basis, miter=1000, tol = 1e-8, rule='bland'):
    m, n = A.shape
    B = LU(A, basis)
    a_indices = np.arange(n)
    b_indices = np.arange(m)

    for iter in range(miter):
        xb = x[basis]
        cb = c[basis]
        bl = np.zeros(n, dtype=bool)
        bl[basis]=1 
        
        v = B.solve(cb, transposed = True)

        ch = c -v.dot(A)
        ch = ch[~bl]
        reduced_cost_pos = np.all(ch >= -tol)
        if reduced_cost_pos:
            break
        pivIndex = select_pivot(ch, bl, a_indices, rule = rule)
        u = B.solve(A[:, pivIndex])
        j = u>tol
        if not np.any(j):
            break
        t = xb[j] / u[j]
        l = np.argmin(t)
        tp = t[l]
        x[basis] = x[basis] - tp * u
        x[pivIndex] = tp
        # Forrest Tomlin Update
        B.update(b_indices[j][l], pivIndex)
        basis = B.b
    return x, basis
# This is Bland's Rule / minimum reduced cost
def select_pivot(ch, bl, a, tol = 1e-8, rule = "bland"):
    rule = rule.lower()
    if rule == "bland":
        return a[~bl][ch < -tol][0]
    elif rule == "minrc":
        return a[~bl][np.argmin(ch)]
    