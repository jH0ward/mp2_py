from numpy import linalg as LA
from numpy import array,zeros
from math import sqrt,pi,exp,erf
import sys

def spoverlap(phi1,phi2,k,p):
 
    AB2=(LA.norm(phi1.center-phi2.center))**2
    alpha1=phi1.alpha[k]
    alpha2=phi2.alpha[p]

    gamma=alpha1+alpha2
    zeta=alpha1*alpha2/float(gamma)
 
     

    Pp=(phi1.center*alpha1+phi2.center*alpha2)/float(gamma)
    PA=Pp-phi1.center
    PB=Pp-phi2.center
    

     # Compute Ix,Iy,and Iz
    l1=phi1.l
    l2=phi2.l

    Ix=fIx(fk,l1,l2,PA[0],PB[0],gamma)

    m1=phi1.m
    m2=phi2.m
 
    Iy=fIx(fk,m1,m2,PA[1],PB[1],gamma)
 
    n1=phi1.n
    n2=phi2.n

    Iz=fIx(fk,n1,n2,PA[2],PB[2],gamma)
 
    result=exp(-1*alpha1*alpha2*AB2/float(gamma))*Ix*Iy*Iz#*\
    return result
 
def factorial(x):
    mult=x-1
    if mult<=0:
        return 1
    result=x
    while mult>0:
        result=result*mult
        mult-=1
    return result

def dfactorial(x):
    if x<=0:
        return 1
    mult=x-2
    if mult<=0:
        return x
    result=x
    while mult>0:
        result=result*mult
        mult=mult-2
    return result


def binomial(n,k):
    if k>n:
        print "Invalid binomial coefficient"
    return factorial(n)/factorial(k)/factorial(n-k)


def fk(k,l1,l2,PAx,PBx):

    if -1*k > k-2*l2:
        q = -1*k
    else:
        q = k-2*l2


    if k < 2*l1-k:
        finish=k
    else:
        finish=2*l1-k

    result=0.0
    start=q
    for l in range(start,finish+1,2):
        i=(k+l)/2
        j=(k-l)/2
        result+=binomial(l1,i)*binomial(l2,j)*PAx**(l1-i)*PBx**(l2-j)

    return result


def fIx(fk,l1,l2,PAx,PBx,gamma):
    result=0.0
    for i in range((l1+l2)/2+1):
        result+=fk(2*i,l1,l2,PAx,PBx)*dfactorial(2*i-1)/float((2*gamma)**i)*\
                    (pi/float(gamma))**(1./2)
    return result

