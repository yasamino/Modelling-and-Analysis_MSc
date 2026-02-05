
def Degree_of_Confluence(f,th,size):
    Area=[]
    m = np.vstack(th.cellInd)
    c = [ f[i][:180] for i in range(len(m)) ]
    for cell in c:
        hull = ConvexHull(cell[:,0:2])
        Area.append(hull.volume)
    return sum(Area)/(size*size)

