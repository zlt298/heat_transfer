import numpy as np
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


class Grid(object):
    def __init__(self,dim_size,mesh_size,k,q_dot):
        """create a grid with the following parameters:
        dim_size: Size of the area as a tuple (L,H)
        mesh_size: matrix dimensions for the mesh (m,n)
        k: Coefficient of conduction
        q_dot: Heat generation
        """
        self.args = (dim_size,mesh_size,k,q_dot)
        self.DIMX,self.DIMY = dim_size
        self.M,self.N = mesh_size
        self.DX,self.DY = self.DIMX/float(self.M),self.DIMY/float(self.N)
        self.K = k
        self.QDOT = q_dot
        
        self.mat_id = [(x,y) for y in range(self.N) for x in range(self.M)]
        self.flat_id = range(len(self.mat_id))
        self.matToflat = dict(zip(self.mat_id,self.flat_id))
        self.matToflat[(-1,-1)] = -1

    def init_mesh(self):
        self.mesh = [Node(x,y,(0,(0,None),0),*self.args) for x,y in self.mat_id]
        self.temperature = np.zeros((self.M,self.N))
        for y in range(self.M):
            for x in range(self.N):
                pass
        
    def __str__(self):
        s = ''
        for x in range(self.M):
            s = s + ' | '.join(['%3.f'%ss for ss in self.temperature[x,:]]) + '\n'
        return s

    def validate_bounds(self):
        return False
    
    def evaluate(self):
        if self.validate_bounds():
            return None
        else:
            print "Wall's are not properly bounded"
            return None

    def show_temp(self):
        return None
    
class Node(Grid):
    def __init__(self,ypos,xpos,bc,*args):
        """create a node with the following parameters:
        xpos: Position in the X, 0 indexed
        ypos: Position in the Y, 0 indexed
        dy: node wall size
        dx: node wall size
        bc:
            The boundary condition for the specified side, given as a tuple:
            (Configuration, Condition, Rotation):
                Configuration
                    internal node = 0, wall node = 1, Internal Corner = 2, External Corner = 3
                Condition
                    None = (0,None)
                    Constant Temperature = (1, T)
                    Constant Flux = (2,q")
                    Convection = (3,h,Tinf)
                Rotation
                    right,up,left,down = 0,1,2,3
        """
        self.DIMX,self.DIMY = args[0]
        self.M,self.N = args[1]
        self.DX,self.DY = self.DIMX/float(self.M),self.DIMY/float(self.N)
        self.K = args[2]
        self.QDOT = args[3]
        
        self.Y,self.X = ypos,xpos
        self.current_temp = 0.
        self.bound = bc
        
        self.update_bc()
        
    def update_bc(self):
        pos = (self.X,self.Y)
        house = [(self.X+1,self.Y),(self.X,self.Y+1),(self.X-1,self.Y),(self.X,self.Y-1)]

        dx,dy,dxsq,dysq = self.DX,self.DY,pow(self.DX,2),pow(self.DY,2)
        rotd = [dysq,dxsq,dysq,dxsq]
        rotdsq = [dy,dx,dy,dx]

        conf = self.bound[0]
        rot = self.bound[2]
        genConst = -self.QDOT*dxsq*dysq/float(self.K)
        posConst = -2*(dxsq+dysq)
        
        #Handle Boundary Conditions, see writeup for derivations
        #Null Condition
        if self.bound[1][0] == 0:
            self.coeff =     [(house[0],            dysq),
                              (house[1],            dxsq),
                              (house[2],            dysq),
                              (house[3],            dxsq),
                              (pos,                 posConst),
                              ((-1,-1),             genConst)]

        #Constant Temperature Condition
        elif self.bound[1][0] == 1:
            self.coeff =     [(self.pos,1),
                              ((-1,-1),             self.bound[1][1])]

        #Constant Flux Condition
        elif self.bound[1][0] == 2:
            q = self.bound[1][1]
            #Wall Configuration
            if conf == 1:
                self.coeff = [(house[rot],          0),
                              (house[(rot+1)%4],    rotdsq[(rot+1)%4]),
                              (house[(rot+2)%4],    2*rotdsq[(rot+2)%4]),
                              (house[(rot+3)%4],    rotdsq[(rot+3)%4]),
                              (pos,                 posConst),
                              ((-1,-1),             genConst-(2*q*rotdsq[rot]/float(self.K)))]
            
            #Internal Corner Configuration
            elif conf == 2:
                self.coeff = [(house[rot],          rotdsq[rot]),
                              (house[(rot+1)%4],    rotdsq[(rot+1)%4]),
                              (house[(rot+2)%4],    2*rotdsq[(rot+2)%4]),
                              (house[(rot+3)%4],    rotdsq[(rot+3)%4]),
                              (pos,                 posConst),
                              ((-1,-1),             genConst-(2*q*rotdsq[rot]/float(self.K)))]
            
            #External Corner Configuration
            elif conf == 3:
                pass

        #Constant Convection Condition
        elif self.bound[1][0] == 3:
            h,Tinf = self.bound[1][1],self.bound[1][2]
            #Wall Configuration
            if conf == 1:
                pass
            
            #Internal Corner Configuration
            elif conf == 2:
                pass
            
            #External Corner Configuration
            elif conf == 3:
                pass

    def set_bc(self,bound):
        self.bound = bound
        self.update_bc()
        
if __name__ == '__main__':
    g = Grid((100,100),(25,25),20,1)
    g.init_mesh()
    for m in range(1000):
        pass
