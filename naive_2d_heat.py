import numpy as np
import matplotlib.pyplot as plt
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
        
        Default bounds are adiabatic walls
        """
        self.args = (dim_size,mesh_size,k,q_dot)
        self.DIMX,self.DIMY = dim_size
        self.M,self.N = mesh_size
        self.DX,self.DY = self.DIMX/float(self.M),self.DIMY/float(self.N)
        self.K = k
        self.QDOT = q_dot
        
        self.mat_id = [(x,y) for x in range(self.M) for y in range(self.N)]
        self.flat_id = range(len(self.mat_id))
        self.matToflat = dict(zip(self.mat_id,self.flat_id))
        self.matToflat[(-1,-1)] = -1

    def init_mesh(self):
        self.mesh = [Node(x,y,(0,(0,None),0),*self.args) for x,y in self.mat_id]
        self.temperature = np.zeros((self.M,self.N))
        for r in self.wall_iter():
            self.mesh[r].set_bc(self.wall_cond(self.mat_id[r],(2,0)))
        
    def __str__(self):
        s = ''
        for y in range(self.N):
            s = s + ' | '.join(['%3.f'%ss for ss in self.temperature[y,:]]) + '\n'
        return s

    def validate_bounds(self):
        return reduce((lambda x,y:x*y),[self.mesh[r].bound[1][0] != 0 for r in self.wall_iter()])
    
    def evaluate(self):
        if self.validate_bounds():
            coeff = np.zeros((self.M*self.N,self.M*self.N))
            answer = np.zeros((self.M*self.N,1))
            for r in self.flat_id:
                temp = np.zeros((1,self.M*self.N))
                for xy,const in self.mesh[r].coeff[:-1]:
                    temp[0,self.matToflat[xy]] = const
                coeff[r,:] = temp
                answer[r] = self.mesh[r].coeff[-1][1]
            self.temperature = np.linalg.solve(coeff,answer)
            #print self.temperature
            self.temperature.shape = (self.M,self.N)
            self.temperature = np.rot90(self.temperature)
            
            return True
        else:
            print "Wall's are not properly bounded"
            return False

    def show_temp(self):
        plt.imshow(self.temperature,interpolation='bilinear', cmap=plt.get_cmap("coolwarm"))
        plt.colorbar()
        plt.show()
        

    def wall_cond(self,pos,condition):
        if   pos == (self.M-1,self.N-1):return (3,condition,0)
        elif pos == (0,self.N-1)       :return (3,condition,1)
        elif pos == (0,0)              :return (3,condition,2)
        elif pos == (self.M-1,0)       :return (3,condition,3)
        else:
            if   pos[0] == self.M-1:return (1,condition,0)
            elif pos[1] == self.N-1:return (1,condition,1)
            elif pos[0] == 0       :return (1,condition,2)
            elif pos[1] == 0       :return (1,condition,3)
            else: return (0,(0,None),0)
    
    def wall_iter(self):
        wall = [self.matToflat[(x,y)] for y in range(self.N) for x in range(self.M) if x == 0 or y == 0 or x == self.M-1 or y == self.N-1]
        for r in wall:
            yield r
    
class Node(Grid):
    def __init__(self,xpos,ypos,bc,*args):
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
        rotd = [dx,dy,dx,dy]
        rotdsq =[dysq,dxsq,dysq,dxsq]

        conf = self.bound[0]
        rot = self.bound[2]
        sadleConst = dx*dy*(dx+dy)
        genConst = -self.QDOT*dxsq*dysq/float(self.K)
        posConst = -2*(dxsq+dysq)
        
        #Handle Boundary Conditions, see writeup for derivations
        #Null Condition=============================================================================
        if self.bound[1][0] == 0:
            self.coeff =     [(house[0],            dysq),
                              (house[1],            dxsq),
                              (house[2],            dysq),
                              (house[3],            dxsq),
                              (pos,                 posConst),
                              ((-1,-1),             genConst)]

        #Constant Temperature Condition=============================================================
        elif self.bound[1][0] == 1:
            self.coeff =     [(pos,                 1),
                              ((-1,-1),             self.bound[1][1])]

        #Constant Flux Condition
        elif self.bound[1][0] == 2:
            q = self.bound[1][1]
            #Wall Configuration ====================================================================
            if conf == 1:
                self.coeff = [(house[rot],          0),
                              (house[(rot+1)%4],    rotdsq[(rot+1)%4]),
                              (house[(rot+2)%4],    2*rotdsq[(rot+2)%4]),
                              (house[(rot+3)%4],    rotdsq[(rot+3)%4]),
                              (pos,                 posConst),
                              ((-1,-1),             genConst-(2*q*rotdsq[rot]*rotd[rot]/float(self.K)))]
            
            #Internal Corner Configuration==========================================================
            elif conf == 2:
                self.coeff = [(house[rot],          rotdsq[rot]),
                              (house[(rot+1)%4],    rotdsq[(rot+1)%4]),
                              (house[(rot+2)%4],    2*rotdsq[(rot+2)%4]),
                              (house[(rot+3)%4],    2*rotdsq[(rot+3)%4]),
                              (pos,                 1.5*posConst),
                              ((-1,-1),             1.5*genConst-(q*(sadleConst)/float(self.K)))]
            
            #External Corner Configuration==========================================================
            elif conf == 3:
                self.coeff = [(house[rot],          0),
                              (house[(rot+1)%4],    0),
                              (house[(rot+2)%4],    rotdsq[(rot+2)%4]),
                              (house[(rot+3)%4],    rotdsq[(rot+3)%4]),
                              (pos,                 0.5*posConst),
                              ((-1,-1),             0.5*genConst-(q*(sadleConst)/float(self.K)))]

        #Constant Convection Condition
        elif self.bound[1][0] == 3:
            h,Tinf = self.bound[1][1],self.bound[1][2]
            #Wall Configuration=====================================================================
            if conf == 1:
                self.coeff = [(house[rot],          0),
                              (house[(rot+1)%4],    rotdsq[(rot+1)%4]),
                              (house[(rot+2)%4],    2*rotdsq[(rot+2)%4]),
                              (house[(rot+3)%4],    rotdsq[(rot+3)%4]),
                              (pos,                 posConst-(2*h*rotdsq[rot]*rotd[rot]/float(self.K))),
                              ((-1,-1),             genConst-(2*h*rotdsq[rot]*rotd[rot]*Tinf/float(self.K)))]
            
            #Internal Corner Configuration==========================================================
            elif conf == 2:
                self.coeff = [(house[rot],          rotdsq[rot]),
                              (house[(rot+1)%4],    rotdsq[(rot+1)%4]),
                              (house[(rot+2)%4],    2*rotdsq[(rot+2)%4]),
                              (house[(rot+3)%4],    2*rotdsq[(rot+3)%4]),
                              (pos,                 1.5*posConst-h*sadleConst/float(self.K)),
                              ((-1,-1),             1.5*genConst-(h*Tinf*(sadleConst)/float(self.K)))]
            
            #External Corner Configuration==========================================================
            elif conf == 3:
                self.coeff = [(house[rot],          0),
                              (house[(rot+1)%4],    0),
                              (house[(rot+2)%4],    rotdsq[(rot+2)%4]),
                              (house[(rot+3)%4],    rotdsq[(rot+3)%4]),
                              (pos,                 0.5*posConst-h*sadleConst/float(self.K)),
                              ((-1,-1),             0.5*genConst-(h*Tinf*(sadleConst)/float(self.K)))]
        self.coeff = [c for c in self.coeff if (c[0]==(-1,-1)) or (c[0][0] != -1 and c[0][0] != self.M and c[0][1] != -1 and c[0][1] != self.N)]
        
    def set_bc(self,bound):
        self.bound = bound
        self.update_bc()
        
if __name__ == '__main__':

    #Initial test
##    g = Grid((50,50),(50,50),20,100)
##    g.init_mesh()
##    
##    for m in range(g.M):
##        g.mesh[m].set_bc((1,(1,300),0))
##        g.mesh[len(g.flat_id) - m - 1].set_bc((1,(1,700),1))
##    g.mesh[g.M-1].set_bc((3,(1,300),0))
##    g.mesh[0].set_bc((3,(1,300),1))
##    g.mesh[len(g.flat_id)-g.M].set_bc((3,(1,700),3))
##    g.mesh[len(g.flat_id)-1].set_bc((3,(1,700),4))
##    
##    g.evaluate()
##    g.show_temp()

##    #Simple wall with heat gen
##    g = Grid((0.2,0.001),(50,50),4,1000)
##    g.init_mesh()
##    
##    for r in g.wall_iter():
##            
##        if g.mat_id[r][0] == g.M-1 and g.mat_id[r][1] != g.N-1 and g.mat_id[r][1] != 0:
##            g.mesh[r].set_bc((1,(3,20,273.15+50),0))
##
##    for r in g.flat_id:
##        if g.mat_id[r][0] < 25 and g.mat_id[r][0] != 0 and g.mat_id[r][1] != 0 and g.mat_id[r][1] != g.N-1:
##            g.mesh[r].QDOT = 0
##            g.mesh[r].update_bc()
##
##    g.evaluate()
##    g.show_temp()

    g = Grid((1,0.01),(51,51),16,900)
    g.init_mesh()

    for r in g.wall_iter():
        if g.mat_id[r][0] == g.M-1:
            g.mesh[r].set_bc(g.wall_cond(g.mat_id[r],(1,307.8540)))#(3,20,296)))
    
##    for r in g.wall_iter():
##        if g.mat_id[r][0] == g.M-1 and g.mat_id[r][1] != g.N-1 and g.mat_id[r][1] != 0:
##            self.mesh[r].set_bc(self.wall_cond(self.mat_id[r],(2,0)))
##            g.mesh[r].set_bc((1,(3,20,293),0))

    for r in g.flat_id:
        if g.mat_id[r][0] < 26:
            g.mesh[r].QDOT = 0
            g.mesh[r].K = 40
            g.mesh[r].update_bc()
    for r in g.flat_id:
        if g.mat_id[r][0] > 39:
            g.mesh[r].QDOT = 0
            g.mesh[r].K = 40
            g.mesh[r].update_bc()

    g.evaluate()
    g.show_temp()
