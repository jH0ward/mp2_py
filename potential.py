from numpy import linalg as LA
from numpy import array,zeros
from math import sqrt,pi,exp,erf
from overlap import factorial
from bvals import *
import sys

def ssAttraction(phi1,phi2,nuc,k,p):
    ''' 
    Function to calculate V integral assoc. with 2 s orbitals
    Usage: value = attraction(phi1,phi2,nuc)
    Input: 2 Gaussian3D instances and 1 list of nuclei coordinates
                ( in separate arrays )
    Output: A float equal to the V integral
    '''

    P=zeros(3)
    C=zeros(3)
    total=0.0

    alpha1=phi1.alpha[k]
    alpha2=phi2.alpha[p]
    gamma=alpha1+alpha2
    zeta=alpha1*alpha2/gamma
    P=(phi1.center*alpha1+phi2.center*alpha2)/gamma
    for i in range(len(nuc)):
        C=array(nuc[i][0])
        PC2=(LA.norm(P-C))**2       
        AB2=(LA.norm(phi1.center-phi2.center))**2
        t=gamma*PC2
        if t==0:  F0=nuc[i][1] # nuc[i][1] holds Z_c for i
        else:
          F0=nuc[i][1]*sqrt(pi)/2/sqrt(t)*erf(sqrt(t))
        result=-2*pi/gamma* \
            exp(-1*zeta*AB2)*F0
        total+=result
    return total
        

def spAttraction(phi1,phi2,nuc,k,p):
    P=zeros(3)
    C=zeros(3)
    total=0.0

    alpha1=phi1.alpha[k]
    alpha2=phi2.alpha[p]
    gamma=alpha1+alpha2
    zeta=alpha1*alpha2/gamma
    P=(phi1.center*alpha1+phi2.center*alpha2)/gamma

    # Determine if this is p_x, p_y, or p_z function
    l1=phi1.l
    l2=phi2.l
    m1=phi1.m
    m2=phi2.m
    n1=phi1.n
    n2=phi2.n

    first='s' # Assume everything is wrong first
    second='s'

    if l1+l2==1:
        comp=0  # Use x-comp.
        if l1==1:
            first='p'
        if l2==1:
            second='p'

    elif m1+m2==1:
        comp=1  # Use y-comp.
        if m1==1:
            first='p'
        if m2==1:
            second='p'
    elif n1+n2==1:
        if n1==1:
            first='p'
        if n2==1:
            second='p'
        comp=2  # Use z-comp.
    else:
        sys.exit()

    for i in range(len(nuc)):
        X_PB=P[comp]-phi2.center[comp] # new
        X_PA=P[comp]-phi1.center[comp] # new
        C=array(nuc[i][0])
        X_PC=P[comp]-C[comp] # new X_PA=P[0]-phi1.center[0] # new PC2=(LA.norm(P-C))**2
        AB2=(LA.norm(phi1.center-phi2.center))**2


        #t=gamma*PC2

        if first=='p':
            total+=X_PA*bigTheta(phi1,phi2,C,0,k,p,nuc,i)\
                    - X_PC*bigTheta(phi1,phi2,C,1,k,p,nuc,i)
                  #+ 1/(2*gamma)*bigTheta(phi1,phi2,C,0,k,p,nuc,i)\
                  #- 1/(2*gamma)*bigTheta(phi1,phi2,C,1,k,p,nuc,i)
        elif second=='p':
            total+=X_PB*bigTheta(phi1,phi2,C,0,k,p,nuc,i)-\
                    X_PC*bigTheta(phi1,phi2,C,1,k,p,nuc,i)#+\
               # 1/(2*gamma)*bigTheta(phi1,phi2,C,0,k,p,nuc,i)-\
               # 1/(2*gamma)*bigTheta(phi1,phi2,C,1,k,p,nuc,i)


    return total


def ppAttraction(phi1,phi2,nuc,k,p):
    P=zeros(3)
    C=zeros(3)
    total=0.0

    alpha1=phi1.alpha[k]
    alpha2=phi2.alpha[p]
    gamma=alpha1+alpha2
    zeta=alpha1*alpha2/gamma
    P=(phi1.center*alpha1+phi2.center*alpha2)/gamma

    l1=phi1.l
    l2=phi2.l
    m1=phi1.m
    m2=phi2.m
    n1=phi1.n
    n2=phi2.n

    if l1+l2==2:
        comp=0
    elif m1+m2==2:
        comp=1
    elif n1+n2==2:
        comp=2
    else:
        return 0.0


    for i in range(len(nuc)):
        X_PB=P[comp]-phi2.center[comp] # new
        X_PA=P[comp]-phi1.center[comp] # new
        C=array(nuc[i][0])
        X_PC=P[comp]-C[comp] # new X_PA=P[0]-phi1.center[0] # new PC2=(LA.norm(P-C))**2
        AB2=(LA.norm(phi1.center-phi2.center))**2


        #t=gamma*PC2

        total+=(X_PB*bigTheta(phi1,phi2,C,0,k,p,nuc,i)-\
                X_PC*bigTheta(phi1,phi2,C,1,k,p,nuc,i))*\
                X_PA\
              + 1/(2*gamma)*bigTheta(phi1,phi2,C,0,k,p,nuc,i)\
              - 1/(2*gamma)*bigTheta(phi1,phi2,C,1,k,p,nuc,i)\
              - (X_PB*bigTheta(phi1,phi2,C,1,k,p,nuc,i)\
              - X_PC*bigTheta(phi1,phi2,C,2,k,p,nuc,i))\
              * X_PC

           # 1/(2*gamma)*bigTheta(phi1,phi2,C,0,k,p,nuc,i)-\
           # 1/(2*gamma)*bigTheta(phi1,phi2,C,1,k,p,nuc,i)

    return total


def coleTAYLOR(x,order,su):
    total=0.0

    for k in range(order+1):
        total+=((-1*x)**k)/factorial(k)/(2*su+2*k+1)
    return total

def mac(x,N,order=6):
    total=0.0
    # Find nearest tenth to x
    xten=x*10
    closest=round(xten)/10.
    diff=x-closest
    for k in range(order+1):
        marker=k+N
        if marker==1:
            mydict=boysONE
        elif marker==2:
            mydict=boysTWO
        elif marker==3:
            mydict=boysTHREE
        elif marker==4:
            mydict=boysFOUR
        elif marker==5:
            mydict=boysFIVE
        elif marker==6:
            mydict=boysSIX
        elif marker==7:
            mydict=boysSEVEN
        elif marker==8:
            mydict=boysEIGHT
        total+=mydict[closest]*(-1*diff)**k/factorial(k)
    return total


def bigTheta(phi1,phi2,C,N,kk,p,nuc,ii):
    alpha1=phi1.alpha[kk]
    Ax=phi1.center[0]
    Ay=phi1.center[1]
    Az=phi1.center[2]

    alpha2=phi2.alpha[p]
    Bx=phi2.center[0]
    By=phi2.center[1]
    Bz=phi2.center[2]

    gamma=alpha1+alpha2
    zeta=alpha1*alpha2/gamma

    const=2*pi/(gamma)
    AB2=(LA.norm(phi1.center-phi2.center))**2
    Kab=exp(-1*zeta*(Ax-Bx)**2)*\
        exp(-1*zeta*(Ay-By)**2)*\
        exp(-1*zeta*(Az-Bz)**2)
    #Kab=exp(-1*zeta*AB2)

    P=(phi1.center*alpha1+phi2.center*alpha2)/gamma


    PC2=(LA.norm(P-C))**2

    t=gamma*PC2
    total=0.0
    if N==0:
        if t<1e-17:
            F=-nuc[ii][1]
        else:
            F=-nuc[ii][1]*sqrt(pi)/2/sqrt(t)*erf(sqrt(t))
        
        total=Kab*const*F

    elif N==1: # N=1
        if t<1e-17:
            F=-nuc[ii][1]*1/float(2*N+1)
        elif t<0.18:
            F=-coleTAYLOR(t,6,N)*nuc[ii][1]
        elif t<19.35:
            F=-nuc[ii][1]*mac(t,N)
        else:
            F=-nuc[ii][1]*0.25*sqrt(pi)/t**1.5

        total=Kab*const*F

    else: #N=2
        if t<1e-17:
            F=-nuc[ii][1]*1/float(2*N+1)
        elif t<0.18:
            F=-coleTAYLOR(t,6,N)*nuc[ii][1]
        elif t<19.35:
            F=-nuc[ii][1]*mac(t,N)
        else:
            F=-nuc[ii][1]*3./8.*sqrt(pi)/t**2.5

        total=Kab*const*F

    return total
