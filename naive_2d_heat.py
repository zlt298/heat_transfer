"""
A simple implementation of a rectangular grid 2D steady state heat transfer 
Grid is sized L x H, and the mesh size is m x n

   X-----------------
Y  [m=1,n=3][m=2,n=3]
|  [m=1,n=2][m=2,n=2]
|  [m=1,n=1][m=2,n=1]

4 boundary conditions on the outside of the grid
Boundary conditions include:
    q = 0
    q = q0
    q = h(Ts-Tinf)
    T = Ts
"""


class grid():
    def __init__(self,dim_size,mesh_size,k,q_dot):
        """create a grid with the following parameters:
        dim_size: Size of the area as a tuple (L,H)
        mesh_size: matrix dimensions for the mesh (m,n)
        k: Coefficient of conduction
        q_dot: Heat generation
        """
        self.DIMX,self.DIMY = dim_size
        self.M,self.N = mesh_size
        self.DX,self.DY = self.DIMX/float(self.M),self.DIMY/float(self.N)
        self.K = k
        self.QDOT = q_dot

    def init_mesh(self):
        self.mesh = [[cell(x,y,(0,(0,None),0)) for x in range(self.M)] for y in range(self.N)]
        self.temperature = np.zeros((self.M,self.N))
        for y in range(self.M):
            for x in range(self.N):
                self.Temperature[x,y] = self.mesh[x][y].current_temp
        
    def __str__(self):
        s = ''
        for x in range(self.M):
            s = s + ' | '.join([str(ss) for ss in self.Temperature[x,:]]) + '\n'
        return s

    def evaluate(self,targeterror = 0.1):
        return None

    def show_temp(self,):
        return None
    
class cell(grid):
    def __init__(self,ypos,xpos,bc):
        """create a cell with the following parameters:
        ypos: Position in the Y, 0 indexed
        xpos: Position in the X, 0 indexed
        dy: cell wall size
        dx: cell wall size
        bc:
            The boundary condition for the specified side, given as a tuple:
            (type,value):
                type:
                    0 = no boundary condition,value as 0 or None
                    1 = constant temperature condition, value as float of the temperature in K
                    2 = constant flux condition, value as float of the flux in watts
                    3 = convection condition, value as tuple, (h,tinf)
            eg. bcTop = (2,-5000),bcRight = (3,(0.6,310)),bcBot = (0,), bcLeft = (0,)
        """
        self.Y,self.X = ypos,xpos
        self.house = [super.mesh[self.Y-1][self.X],super.mesh[self.Y][self.X-1],
                  super.mesh[self.Y+1][self.X],super.mesh[self.Y][self.X+1]]
        
        self.current_temp = 0.
        self.prev_temp = 0.
        self.error = 1e10
        self.bound = (bcTop,bcLeft,bcBot,bcRight)
        self.bounds = sum(val[0]>0 for val in self.bound)
        self.updateBC()
        
    def updateBC():
        self.CONSTANTTEMP = sum([t[1] for t in self.bound if t[0]==1])
        if self.CONSTANTTEMP:self.CONSTANTTEMP = self.CONSTANTTEMP/float(len([None for t in self.bound if t[0]==1]))

        dxsq,dysq = pow(super.DX,2),pow(super.DY,2)
        self.equation = lambda z:(dysq*(self.house[1].current_temp+self.house[3].current_temp)+
                                  dxsq*(self.house[0].current_temp+self.house[2].current_temp)+
                                  super.QDOT/float(super.K)*dysq*dxsq)/float(2*(dxsq+dysq))
        #Handle Boundary Conditions
        if self.CONSTANTTEMP: #Temperature Constraint overrides all else
            self.current_temp = self.CONSTANTTEMP
        elif self.bounds == 1:
            dxy = [(super.DY,super.DX),(super.DX,super.DY),(super.DY,super.DX),(super.DX,super.DY)]
            pos = [val[0]>0 for x in a].index(1)
            if sum(val[0] for val in self.bound) == 2:   # Plane Wall q
                self.equation = lambda z:(super.DY*(self.neighbor[1]+self.neighbor[3])+super.DX*(self.neighbor[0]+self.neighbor[2])+
                                          super.QDOT/float(super.K)*pow(super*DY,2)*pow(super.DX,2))/float(2*(super.DX+super.DY))
            else:                                        # Plane Wall h
                pass
        elif self.bounds == 2:
            if sum(val[0] for val in self.bound) == 4:   # Corner q-q
                pass
            elif sum(val[0] for val in self.bound) == 6: # Corner h-h
                pass
            else:                                        # Corner h-q
                pass
        elif self.bounds == 0:
            pass
        else:
            print 'Cell is overconstrained. Check (%s,%s)'%(self.Y,self.X)

    def setBC(side,val):
        tempbound = list(self.bound)
        tempbound[side] = val
        self.bound = tuple(tempbound)
        self.bounds = sum(val[0]>0 for val in self.bound)
        self.updateBC()
        
if __name__ == '__main__':
    g = grid((100,100),(1000,1000),20,0)
    for m in range(1000):
        pass
