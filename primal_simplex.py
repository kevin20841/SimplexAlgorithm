import numpy as np
import sys
np.set_printoptions(threshold=sys.maxsize)

# First Implementation: Primal Simplex, with tableaus. This should be the most inefficient.
    # Phase 1 objective
     
def init_tableau_ID(A, b, c):
    [n,m] = A.shape
    n_contraint = np.less(b, 0)
    A[n_contraint] *=-1
    b[n_contraint] *=-1
    
    row_con = np.hstack((A, np.eye(n), b[:, np.newaxis]))
    av = np.arange(n) + m
    
    basis = av.copy()
    row_obj = np.hstack((c, np.zeros(n+1)))
    row_pseudo_obj = - row_con.sum(axis = 0)
    row_pseudo_obj[av] = 0
    tab = np.vstack((row_con, row_obj, row_pseudo_obj))
    
    return tab, basis, av

def can_improve(tab):
    return np.any(tab[-1, :-1]>0)

# Bland's Rule at first. Choose first non-basic variable as entering variable
# np.where(t>0, t, np.inf).argmin()
def choose_pivot_col(tab, tol = 1e-10, bland=False):
    ma = np.ma.masked_where(tab[-1, :-1] >= -tol, tab[-1, :-1], copy=False)
    if ma.count() == 0:
        return False, np.nan
    if bland:
        return True, np.where(tab[-1] > 0)[0][0]
    return True, np.ma.nonzero(ma == ma.min())[0][0]


# This is pricing! 
def choose_pivot_row(tab, pivColi, basis, phase, tol = 1e-10, bland=False):
    if phase == 1:
        k = 2
    else:
        k = 1
    
    ma = np.ma.masked_where(tab[:-k, pivColi] <= tol, tab[:-k, pivColi], copy=False)
    if ma.count() == 0:
        return False, np.nan
    mb = np.ma.masked_where(tab[:-k, pivColi] <= tol, tab[:-k, -1], copy=False)

    residuals = mb/ma
    min_rows = np.ma.nonzero(residuals == residuals.min())[0]
    if bland:
        return True, min_rows[np.argmin(np.take(basis, min_rows))]
    return True, min_rows[0]

def apply_pivot(tab, basis, pivRowi, pivColi):
    basis[pivRowi] = pivColi
    tab[pivRowi, :] = tab[pivRowi] / tab[pivRowi][pivColi]
    for i in range(tab.shape[0]):
        if i != pivRowi:
            tab[i] = tab[i] - tab[pivRowi]  * tab[i, pivColi]
    
# Requires two phases
# Add slack variables, then find a starting feasible solution by adding a pseudo-objective
# Remove the pseudo-objective and solve the problem for real. 

def primal_simplex(A, b, c, tol=1e-8):
    A2 = np.eye(A.shape[0])

    c = np.concatenate((c, np.zeros((A.shape[0],))))
    A = np.hstack([A, A2])
    tab, basis, av = init_tableau_ID(A, b, c)
    # Phase 1:
    solve_simplex(tab, basis, phase = 1, tol=tol)
    if (abs(tab [-1, -1]) < tol):
        tab = tab[:-1, :]
        tab = np.delete(tab, av, 1)
    else:
        print("There was no feasible solution found in phase 1.")
        return None, None
    #Found starting point, doing phase 2
    solve_simplex(tab, basis, phase=2, tol=tol)
    n, m= A.shape
    solution = np.zeros(m+n)
    solution[basis[:n]] = tab[:n, -1]
    return solution[:m], tab[-1, -1]

# Pseudocode for simplex is:
# 1. Generate the Initial Tableau
# 2. While the tableau can be improved:
# 3.    Find the pivot index
# 4.    Pivot about the tableau
# 5. Return both info from tableau and the optimal value

def solve_simplex(tab, basis, phase, tol = 10e-10):
    if phase == 2:
        for pivrow in [row for row in range(basis.size)
                       if basis[row] > tab.shape[1] - 2]:
            non_zero_row = [col for col in range(tab.shape[1] - 1)
                            if abs(tab[pivrow, col]) > tol]
            if len(non_zero_row) > 0:
                pivcol = non_zero_row[0]
                apply_pivot(tab, basis, pivrow, pivcol)
    
    finished = False
    bland = False
    count = 0
    while (not finished):
        c1, pivColi = choose_pivot_col(tab, bland=bland, tol=tol)
        if not c1:
            break
        c2, pivRowi = choose_pivot_row(tab, pivColi, basis,phase, bland=bland, tol=tol)
        if not c2:
            break

        apply_pivot(tab, basis, pivRowi, pivColi)

        count+=1

