"""
A simple implementation of a rectangular grid 2D steady state heat transfer 
Grid is sized L x H, and the mesh size is m x n

   X-----------------
Y  [m=1,n=1][m=2,n=1]
|  [m=1,n=2][m=2,n=2]
|  [m=1,n=3][m=2,n=3]

4 boundary conditions on the outside of the grid
Boundary conditions include:
    q = 0
    q = q0
    q = h(Ts-Tinf)
    T = Ts
"""
import numpy as np
import matplotlib as mpl

class grid(object,dim_size,mesh_size,k,q_dot):

    def __init__(self):
        pass

    def __str__(self):
        pass

class point(grid):
    def __init__(self):
        pass

    def __str__(self):
        pass
