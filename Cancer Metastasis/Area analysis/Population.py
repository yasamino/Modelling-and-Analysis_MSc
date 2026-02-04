import numpy as np

def count_cells(f,th):
    m = np.vstack(th.cellInd)
    total = len(m)
    s1 = np.where(m < 0)[0]
    s2 = np.where(m >= 0 )[0]
    h = [ f[i][:180] for i in s1 ]
    soft = len(h)  
    g = [ f[i][:180] for i in s2 ]
    hard = len(g)
    return soft, hard, total