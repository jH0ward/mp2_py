import numpy as np
from math import *
import sys


def new(bf):
    nprims=bf.nprims
    l=bf.l
    m=bf.m
    n=bf.n
    total=0.0
    cc=bf.cc
    L=l+m+n
    alpha=bf.alpha
    for i in range(nprims):
        for j in range(nprims):
            #total+=cc[i]*cc[j]*(2*alpha[i]/pi)**0.75*(2*alpha[j]/pi)**0.75 \
            #/float((alpha[i]+alpha[j])**(1.5))
            total+=cc[i]*cc[j]*((2/pi)**1.5)*2**(l+m+n)\
                    *alpha[i]**((2*l+2*m+2*n+3)/4.0)\
                    *alpha[j]**((2*l+2*m+2*n+3)/4.0)\
                    /dfact(2*l-1)/dfact(2*m-1)/dfact(2*n-1)\
            /float((alpha[i]+alpha[j])**(L+1.5))

     

    
    #N =(np.pi**1.5 * dfact(2*l-1)*dfact(2*m-1)*dfact(2*n-1) / \
    #         2**(l+m+n)*total)**(-0.5)

    #N=1.0/sqrt(np.pi**1.5*total)
    norm=total*pi**(1.5)
    cnorm=1.0/sqrt(norm)


    bf.NC=cnorm
    return 0



def dfact(x):
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
