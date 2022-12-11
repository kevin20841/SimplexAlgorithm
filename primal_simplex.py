import numpy as np
np.set_printoptions(formatter={'float': lambda x: "{0:0.3f}".format(x)})
# First Implementation: Primal Simplex, with tableaus. This should be the most inefficient.
def init_tableau_ID(A, b, c):
    [m,n] = A.shape
    ret = np.zeros((m+1, m+n+2))
    ret[:m, :n] = A
    ret[:m+1, n:n+m+1] = np.identity(m+1)
    ret[m, :n] = -c
    ret[:m, -1] = b
    return ret

def can_improve(tab):
    return not (np.all(tab[-1] >=0))

# Bland's Rule at first. Choose first non-basic variable as entering variable
# np.where(t>0, t, np.inf).argmin()
def choose_pivot_col(tab):
    return np.where(tab[-1] != 0)[0][0]
    #return np.argmin(tab[-1])

# This is pricing! 
def choose_pivot_row(tab, pivColi):
    pivCol = tab[:, pivColi]


    residuals = np.divide(tab[:, -1], pivCol, out = np.full(pivCol.shape, np.inf), where = (pivCol > 0))

    if np.all(pivCol <=0):
        return "Unbounded Minimum"
    return residuals.argmin()

    
# General Pseudocode is:
# 1. Generate the Initial Tableau
# 2. While the tableau can be improved:
# 3.    Find the pivot index
# 4.    Pivot about the tableau
# 5. Return both info from tableau and the optimal value
# 
# Indices:
# obj = -tab[m, n]
# tb = tab[:m, n]
# tc = tab[m, :n]
# tA = tab[:m, :n]

def primal_simplex(A, b, c):
    [m, n] = A.shape
    tab = init_tableau_ID(A, b, c)
    temp = 0
    while (can_improve(tab)):
        pivColi = choose_pivot_col(tab)
        pivRowi = choose_pivot_row(tab, pivColi)
        if pivRowi =="Unbounded Minimum":
            return "Unbounded Minimum"
        tab[:, pivColi] = tab[:, pivColi] / tab[pivRowi][pivColi]
        i = 0
        while i < m+1:
            if i != pivColi:
                tab[i] = tab[i] - tab[pivRowi] * tab[i][pivColi]
            i+=1
    # TODO Solve for the feasible solution
    return tab, tab[m, m+n+1]
