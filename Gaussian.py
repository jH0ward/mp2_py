from numpy import array
from math import sqrt,pi

class Gaussian3D:
    '''
    Class of 3-d Gaussian-type orbitals.
    Methods:
        .normalize() assigns proper normalization constant to orbital
        .scaleCGS(zeta) scales the default exponent and contraction
                coefficient based on user input (zeta)
    Attributes:
        self.cc - contraction coefficient
        self.Ax - x coordinate of Gaussian orbital's center
        self.Ay - y coordinate of Gaussian orbital's center
        self.Az - z coordinate of Gaussian orbital's center
        self.NC - normalization constant of Gaussian orbital (not used )
        self.center - array of [Ax,Ay,Az]
        self.label - distinguishes Gaussian instance with string (if desired)
    History:
        Updated on 2/29/2012
    '''
    #def __init__(self,symbol,coords,cc,alpha,label):
    def __init__(self,symbol,coords,basis,am,label,mark):
        self.Ax=coords[0]
        self.Ay=coords[1]
        self.Az=coords[2]
        #self.NC=1.0
        self.symbol=symbol
        self.center=array([self.Ax,self.Ay,self.Az])
        self.coords=coords
        self.label=label
        self.basis=basis
        self.mark=mark


        # New code
        self.am=am
        self.nprims=len(basis[1])

        self.l=0
        self.m=0
        self.n=0

        if am[0]==1:
            self.l=1

        if am[1]==1:
            self.m=1

        if am[2]==1:
            self.n=1

        self.L=self.l+self.m+self.n
        self.cc=[]
        self.alpha=[]
        for i in range(self.nprims):
            self.alpha.append(basis[mark][i][0])
            self.cc.append(basis[mark][i][1])

