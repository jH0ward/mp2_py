from Gaussian import Gaussian3D
from overlap import spoverlap
from numpy import linalg as LA
from numpy import array,zeros
from math import sqrt,pi,exp,erf
import sys

def akinetic(phi1,phi2,k,p):
    '''
    Ix = 1/2 l1*l2 <-1|-1>x + 2*alpha1*alpha2* <+1|+1>x
        -alpha1*l2 * <+1|-1>x - alpha2*l1*<-1|+1>x
    '''
    nprims1=phi1.nprims
    nprims2=phi2.nprims
    alpha1=phi1.alpha[k]
    alpha2=phi2.alpha[p]
    # Do term 1 unless one of the basis functions is 'S'

    if phi1.l==0 or phi2.l==0:
        term1=0.0 
    else:
    # Create new gaussian instances with l lowered by 1
        tmpGau1=Gaussian3D(phi1.symbol,phi1.coords,\
                    phi1.basis,phi1.am,-1*phi1.label,phi1.mark)
        tmpGau1.l=phi1.l-1

        tmpGau2=Gaussian3D(phi2.symbol,phi2.coords,\
                    phi2.basis,phi2.am,-1*phi2.label,phi2.mark)
        tmpGau2.l=phi2.l-1

        tmpGau1.nprims=1
        tmpGau2.nprims=1
        tmpGau1.alpha[0]=alpha1
        tmpGau2.alpha[0]=alpha2
    # Comput their overlap and define term1
        term1=1./2*phi1.l*phi2.l*\
             spoverlap(tmpGau1,tmpGau2,k,p)

    # Do term 2 no matter what 
    tmpGau1=Gaussian3D(phi1.symbol,phi1.coords,\
             phi1.basis,phi1.am,-1*phi1.label,phi1.mark)

    tmpGau2=Gaussian3D(phi2.symbol,phi2.coords,\
             phi2.basis,phi2.am,-1*phi2.label,phi2.mark)

    tmpGau1.l=phi1.l+1
    tmpGau2.l=phi2.l+1

    tmpGau1.nprims=1
    tmpGau2.nprims=1

    tmpGau1.alpha[0]=alpha1
    tmpGau2.alpha[0]=alpha2



    term2=2*alpha1*alpha2*\
            spoverlap(tmpGau1,tmpGau2,k,p)
    # Do terunless phi2 is 'S'
    if phi2.l==0:
        term3=0.0
    else:
        tmpGau1.l=phi1.l+1
        tmpGau2.l=phi2.l-1
    
        term3=-1*alpha1*phi2.l*\
                spoverlap(tmpGau1,tmpGau2,k,p)

    # Do terunless phi1 is 'S'
    if phi1.l==0:
        term4=0.0
    else:
        tmpGau1.l=phi1.l-1
        tmpGau2.l=phi2.l+1

        term4=-1*alpha2*phi1.l*\
                spoverlap(tmpGau1,tmpGau2,k,p)
    
    Ix=(term1+term2+term3+term4)


####  Begin m angular momentum

    if phi1.m==0 or phi2.m==0:
        term1=0.0 
    else:
    # Create new gaussian instances with l lowered by 1
        tmpGau1=Gaussian3D(phi1.symbol,phi1.coords,\
                    phi1.basis,phi1.am,-1*phi1.label,phi1.mark)
        tmpGau1.m=phi1.m-1

        tmpGau2=Gaussian3D(phi2.symbol,phi2.coords,\
                    phi2.basis,phi2.am,-1*phi2.label,phi2.mark)
        tmpGau2.m=phi2.m-1

        tmpGau1.nprims=1
        tmpGau2.nprims=1
        tmpGau1.alpha[0]=alpha1
        tmpGau2.alpha[0]=alpha2
    # Computeir overlap and define term1
        term1=1./2*phi1.m*phi2.m*\
             spoverlap(tmpGau1,tmpGau2,k,p)

    # Do term 2 no matter what (reuse tmpGau1 and 2)
    tmpGau1=Gaussian3D(phi1.symbol,phi1.coords,\
              phi1.basis,phi1.am,-1*phi1.label,phi1.mark)

    tmpGau2=Gaussian3D(phi2.symbol,phi2.coords,\
              phi2.basis,phi2.am,-1*phi2.label,phi2.mark)

    tmpGau1.m=phi1.m+1
    tmpGau2.m=phi2.m+1

    tmpGau1.nprims=1
    tmpGau2.nprims=1

    tmpGau1.alpha[0]=alpha1
    tmpGau2.alpha[0]=alpha2


    term2=2*alpha1*alpha2*\
            spoverlap(tmpGau1,tmpGau2,k,p)

    # Do the third unless phi2 is 'S'
    if phi2.m==0:
        term3=0.0
    else:
        tmpGau1.m=phi1.m+1
        tmpGau2.m=phi2.m-1
    
        term3=-1*alpha1*phi2.m*\
                spoverlap(tmpGau1,tmpGau2,k,p)

    # Do term 4 unless phi1 is 'S'
    if phi1.m==0:
        term4=0.0
    else:
        tmpGau1.m=phi1.m-1
        tmpGau2.m=phi2.m+1

        term4=-1*alpha2*phi1.m*\
                spoverlap(tmpGau1,tmpGau2,k,p)
    
    Iy=(term1+term2+term3+term4)


####### Begin n angular momentum

    if phi1.n==0 or phi2.n==0:
        term1=0.0 
    else:
    # Create new gaussian instances with l lowered by 1
        tmpGau1=Gaussian3D(phi1.symbol,phi1.coords,\
                    phi1.basis,phi1.am,-1*phi1.label,phi1.mark)
        tmpGau1.n=phi1.n-1

        tmpGau2=Gaussian3D(phi2.symbol,phi2.coords,\
                    phi2.basis,phi2.am,-1*phi2.label,phi2.mark)
        tmpGau2.n=phi2.n-1

        tmpGau1.nprims=1
        tmpGau2.nprims=1
        tmpGau1.alpha[0]=alpha1
        tmpGau2.alpha[0]=alpha2

    # Comput their overlap and define term1
        term1=1./2*phi1.n*phi2.n*\
             spoverlap(tmpGau1,tmpGau2,k,p)

    # Do term 2 no matter what (reuse tmpGau1 and 2)
    tmpGau1=Gaussian3D(phi1.symbol,phi1.coords,\
              phi1.basis,phi1.am,-1*phi1.label,phi1.mark)

    tmpGau2=Gaussian3D(phi2.symbol,phi2.coords,\
              phi2.basis,phi2.am,-1*phi2.label,phi2.mark)

    tmpGau1.n=phi1.n+1
    tmpGau2.n=phi2.n+1

    tmpGau1.nprims=1
    tmpGau2.nprims=1

    tmpGau1.alpha[0]=alpha1
    tmpGau2.alpha[0]=alpha2


    term2=2*alpha1*alpha2*\
            spoverlap(tmpGau1,tmpGau2,k,p)


    # Do term 3 unless phi2 is 'S'
    if phi2.n==0:
        term3=0.0
    else:
        tmpGau1.n=phi1.n+1
        tmpGau2.n=phi2.n-1
    
        term3=-1*alpha1*phi2.n*\
                spoverlap(tmpGau1,tmpGau2,k,p)

    # Do term 4 unless phi1 is 'S'
    if phi1.n==0:
        term4=0.0
    else:
        tmpGau1.n=phi1.n-1
        tmpGau2.n=phi2.n+1

        term4=-1*alpha2*phi1.n*\
                spoverlap(tmpGau1,tmpGau2,k,p)
    Iz=(term1+term2+term3+term4)
    result=Ix+Iy+Iz

    return result
