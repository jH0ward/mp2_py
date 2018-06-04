from numpy import linalg as LA
from numpy import array,zeros
from math import sqrt,pi,exp,erf
from overlap import factorial
from bvals import *
import sys

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
        if marker>7:
            return total
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



def ssssRepulsion(pair1,pair2,myprims):
    '''
    Function to calculate two-electron repulsion integral for 4 s orbitals
    Usage: value = repulsion(phis)
    Input: 1 list of 4 Gaussian3D instances
    Output: A float equal to the value of the J (or K) integral
    '''
    phi1=pair1[0]
    phi2=pair1[1]
    phi3=pair2[0]
    phi4=pair2[1]

    p=myprims[0]
    q=myprims[1]
    r=myprims[2]
    s=myprims[3]

    result=0.0

    alpha1=phi1.alpha[p]
    alpha2=phi2.alpha[q]
    alpha3=phi3.alpha[r]
    alpha4=phi4.alpha[s]

    cc1=phi1.cc[p]
    cc2=phi2.cc[q]
    cc3=phi3.cc[r]
    cc4=phi4.cc[s]
    gammaP=alpha1+alpha2
    gammaQ=alpha3+alpha4

    zetaP=alpha1*alpha2/gammaP
    zetaQ=alpha3*alpha4/gammaQ

    xi=gammaP*gammaQ/(gammaP+gammaQ)*1.0
    
    P=(alpha1*phi1.center+alpha2*phi2.center)/gammaP
    Q=(alpha3*phi3.center+alpha4*phi4.center)/gammaQ
    
    AB2=(LA.norm(phi1.center-phi2.center))**2
    CD2=(LA.norm(phi3.center-phi4.center))**2
    PQ2=(LA.norm(P-Q))**2

    Kab=exp(-1*zetaP*AB2)
    Kcd=exp(-1*zetaQ*CD2)


    t=xi*PQ2
    if t==0:
        F0=1.0
    else:
        F0=sqrt(pi)/2/sqrt(t)*erf(sqrt(t))
    result+=2*pi**(5./2.)*\
            Kab*Kcd/gammaP/gammaQ/sqrt(gammaP+gammaQ)*F0#*\
        #(2*alpha1/pi)**(3./4.)*(2*alpha2/pi)**(3./4.)*\
        #(2*alpha3/pi)**(3./4.)*(2*alpha4/pi)**(3./4.)
                                                         
                    

    return result

def ssspRepulsion(pair1,pair2,myprims):
    tphi1=pair1[0]
    tphi2=pair1[1]
    tphi3=pair2[0]
    tphi4=pair2[1]

    tp=myprims[0]
    tq=myprims[1]
    tr=myprims[2]
    ts=myprims[3]

    # Reorder phis
    if tphi4.L>0:
        phi1=tphi1
        p=tp
        phi2=tphi2
        q=tq
        phi3=tphi3
        r=tr
        phi4=tphi4
        s=ts
        # do nothin
    elif tphi3.L>0:
        phi3=tphi4
        r=ts
        phi4=tphi3
        s=tr
        phi1=tphi1
        p=tp
        phi2=tphi2
        q=tq

    elif tphi2.L>0:
        phi2=tphi4
        q=ts
        phi4=tphi2
        s=tq
        phi1=tphi3
        p=tr
        phi3=tphi1
        r=tp

    elif tphi1.L>0:
        phi1=tphi4
        p=ts
        phi4=tphi1
        s=tp
        phi2=tphi3
        q=tr
        phi3=tphi2
        r=tq


    alpha1=phi1.alpha[p]
    alpha2=phi2.alpha[q]
    alpha3=phi3.alpha[r]
    alpha4=phi4.alpha[s]


    gammaP=alpha1+alpha2
    gammaQ=alpha3+alpha4
    xi=gammaP*gammaQ/(gammaP+gammaQ)


    A=phi1.center
    B=phi2.center
    C=phi3.center
    D=phi4.center


    P=(alpha1*A+alpha2*B)/gammaP
    Q=(alpha3*C+alpha4*D)/gammaQ
    Q=(C*alpha3+D*alpha4)/gammaQ


    ## Find out the orienation of the one p function
    lxTOT=phi1.l+phi2.l+phi3.l+phi4.l
    lyTOT=phi1.m+phi2.m+phi3.m+phi4.m
    lzTOT=phi1.n+phi2.n+phi3.n+phi4.n

    ## In next block,mycomp is vector component
    ##    and wc is 'which center' (A,B,C,D)


    if phi4.l>0:
        mycomp=0
    elif phi4.m>0:
        mycomp=1
    elif phi4.n>0:
        mycomp=2
    else:
        return 0.0


    XPD=P[mycomp]-D[mycomp]
    XQP=Q[mycomp]-P[mycomp]
    XPQ=P[mycomp]-Q[mycomp]
    XQD=Q[mycomp]-D[mycomp]

   
    term1=XQD*terTheta(phi1,phi2,phi3,phi4,0,p,q,r,s)
    term2=-1*xi*XQP/gammaQ*terTheta(phi1,phi2,phi3,phi4,1,p,q,r,s)

    return term1+term2

    
def ssppRepulsion(pair1,pair2,myprims):
    tphi1=pair1[0]
    tphi2=pair1[1]
    tphi3=pair2[0]
    tphi4=pair2[1]

    tp=myprims[0]
    tq=myprims[1]
    tr=myprims[2]
    ts=myprims[3]

    # Reorder phis
    if tphi4.L>0 and tphi3.L>0:
        phi1=tphi3
        p=tr
        phi2=tphi4
        q=ts
        phi3=tphi1
        r=tp
        phi4=tphi2
        s=tq

    elif tphi1.L>0 and tphi2.L>0:
        phi1=tphi1
        p=tp
        phi2=tphi2
        q=tq
        phi3=tphi3
        r=tr
        phi4=tphi4
        s=ts

    alpha1=phi1.alpha[p]
    alpha2=phi2.alpha[q]
    alpha3=phi3.alpha[r]
    alpha4=phi4.alpha[s]


    gammaP=alpha1+alpha2
    gammaQ=alpha3+alpha4
    xi=gammaP*gammaQ/(gammaP+gammaQ)


    A=phi1.center
    B=phi2.center
    C=phi3.center
    D=phi4.center

    P=(alpha1*A+alpha2*B)/gammaP
    Q=(alpha3*C+alpha4*D)/gammaQ

    # Determine if phi3 and phi4 are x,y,or z

    if phi1.l>0:
        mycomp1=0
    elif phi1.m>0:
        mycomp1=1
    elif phi1.n>0:
        mycomp1=2

    if phi2.l>0:
        mycomp2=0
    elif phi2.m>0:
        mycomp2=1
    elif phi2.n>0:
        mycomp2=2
    else:
        return 0.0

    XPC=P[mycomp1]-C[mycomp1]
    XPQ=P[mycomp1]-Q[mycomp1]

    XPD=P[mycomp2]-D[mycomp2]
    XPQ2=P[mycomp2]-Q[mycomp2]

    XPA=P[mycomp1]-A[mycomp1]
    XAB=A[mycomp1]-B[mycomp1]
    T00=terTheta(phi1,phi2,phi3,phi4,0,p,q,r,s)
    T01=terTheta(phi1,phi2,phi3,phi4,1,p,q,r,s)
    T02=terTheta(phi1,phi2,phi3,phi4,2,p,q,r,s)

    T1000=XPA*T00-xi/gammaP*XPQ*T01
    T1000plus=XPA*T01-xi/gammaP*XPQ*T02
    T2000=XPA*T1000-xi/gammaP*XPQ*T1000plus\
        +0.5/gammaP*(T00-xi/gammaP*T01)
    T1100=T2000+XAB*T1000
    if phi1.am!=phi2.am:
        YPB=P[mycomp2]-B[mycomp2]
        YPQ=P[mycomp2]-Q[mycomp2]
        return YPB*T1000-xi/gammaP*YPQ*T1000plus

    result=T1100


    return result

def spspRepulsion(pair1,pair2,myprims):
    tphi1=pair1[0]
    tphi2=pair1[1]
    tphi3=pair2[0]
    tphi4=pair2[1]

    tp=myprims[0]
    tq=myprims[1]
    tr=myprims[2]
    ts=myprims[3]

    #working copy #if tphi4.L>0 and tphi2.L>0:
        #phi4=tphi4
        #s=ts
        #phi2=tphi2
        #q=tq
        #phi3=tphi3
        #r=tr
        #phi1=tphi1
        #p=tp
    if tphi1.L>0 and tphi3.L>0:
        phi1=tphi1
        p=tp
        phi2=tphi2
        q=tq
        phi3=tphi3
        r=tr
        phi4=tphi4
        s=ts

    elif tphi1.L>0 and tphi4.L>0:
        phi1=tphi1
        p=tp
        phi2=tphi2
        q=tq
        phi3=tphi4
        r=ts
        phi4=tphi3
        s=tr

    elif tphi2.L>0 and tphi3.L>0:
        phi1=tphi2
        p=tq
        phi2=tphi1
        q=tp
        phi3=tphi3
        r=tr
        phi4=tphi4
        s=ts

    elif tphi2.L>0 and tphi4.L>0:
        phi1=tphi2
        p=tq
        phi2=tphi1
        q=tp
        phi3=tphi4
        r=ts
        phi4=tphi3
        s=tr

   # elif tphi4.L>0 and tphi3.L>0:
   #     phi4=tphi4
   #     s=ts
   #     phi2=tphi3
   #     q=tr
   #     phi3=tphi2
   #     r=tq
   #     phi1=tphi1
   #     p=tp
   # elif tphi4.L>0 and tphi1.L>0:
   #     phi4=tphi4
   #     s=ts
   #     phi2=tphi1
   #     q=tp
   #     phi3=tphi2
   #     r=tq
   #     phi1=tphi3
   #     p=tr
   # elif tphi3.L>0 and tphi2.L>0:
   #     phi4=tphi3
   #     s=tr
   #     phi3=tphi4
   #     r=ts
   #     phi2=tphi2
   #     q=tq
   ##     phi1=tphi1
   #     p=tp
   # elif tphi3.L>0 and tphi1.L>0:
   #     phi4=tphi3
   #     s=tr
   #     phi3=tphi4
   #     r=ts
   #     phi2=tphi1
   #     q=tp
   #     phi1=tphi2
   #     p=tq
    #elif tphi2.L>0 and tphi1.L>0:
    #    phi2=tphi2
    #    q=tq
    #    phi4=tphi1
    #    s=tp
    #    phi3=tphi3
    #    r=tr
    #    phi1=tphi4
    #    p=ts
    # Reorder phis
    #if tphi4.L>0:
    #    phi1=tphi1
    #    p=tp
    #    phi2=tphi2
    #    q=tq
    #    phi3=tphi3
    #    r=tr
    #    phi4=tphi4
    #    s=ts
    #    # do nothin
    #    pass
    #elif tphi3.L>0:
    #    phi3=tphi4
    #    r=ts
    #    phi4=tphi3
    #    s=tr
    #    phi1=tphi1
    #    p=tp
    #    phi2=tphi2
    #    q=tq
#
#    elif tphi2.L>0:
#        phi2=tphi4
#        q=ts
#        phi4=tphi2
#        s=tq
#        phi1=tphi1
#        p=tp
#        phi3=tphi3
#        r=tr
#
#    elif tphi1.L>0:
#        phi1=tphi4
#        p=ts
#        phi4=tphi1
#        s=tp
#        phi2=tphi2
#        q=tq
#        phi3=tphi3
#        r=tr
#    #phi1=pair1[0]
#    #phi2=pair1[1]
#    #phi3=pair2[0]
#    #phi4=pair2[1]
##
##    p=myprims[0]
##    q=myprims[1]
##    r=myprims[2]
#    s=myprims[3]

    alpha1=phi1.alpha[p]
    alpha2=phi2.alpha[q]
    alpha3=phi3.alpha[r]
    alpha4=phi4.alpha[s]


    gammaP=alpha1+alpha2
    gammaQ=alpha3+alpha4
    xi=gammaP*gammaQ/(gammaP+gammaQ)


    A=phi1.center
    B=phi2.center
    C=phi3.center
    D=phi4.center

    P=(alpha1*A+alpha2*B)/gammaP
    Q=(alpha3*C+alpha4*D)/gammaQ

    # Determine if phi2 and phi4 are x,y,or z
    if phi1.l>0:
        mycomp1=0
    elif phi1.m>0:
        mycomp1=1
    elif phi1.n>0:
        mycomp1=2

    if phi3.l>0:
        mycomp2=0
    elif phi3.m>0:
        mycomp2=1
    elif phi3.n>0:
        mycomp2=2

    XPA=P[mycomp1]-A[mycomp1]
    XPQ=P[mycomp1]-Q[mycomp1]
    YPQ=P[mycomp2]-Q[mycomp2]
    YQC=Q[mycomp2]-C[mycomp2]
    T00=terTheta(phi1,phi2,phi3,phi4,0,p,q,r,s)
    T01=terTheta(phi1,phi2,phi3,phi4,1,p,q,r,s)
    T02=terTheta(phi1,phi2,phi3,phi4,2,p,q,r,s)

    T1000=XPA*T00-xi/gammaP*XPQ*T01
    T1000plus=XPA*T01-xi/gammaP*XPQ*T02
    if mycomp1!=mycomp2:
        return YQC*T1000+xi/gammaQ*YPQ*T1000plus

    XAB=A[mycomp1]-B[mycomp1]
    XCD=C[mycomp1]-D[mycomp1]
    XPQ=P[mycomp1]-Q[mycomp1]
    XPA=P[mycomp1]-A[mycomp1]

    T1000plus=XPA*T01-xi/gammaP*XPQ*T02
    T2000=XPA*T1000-xi/gammaP*XPQ*T1000plus\
       +0.5/gammaP*(T00-xi/gammaP*T01)

    result=-1*(alpha2*XAB+alpha4*XCD)/gammaQ*T1000\
        +0.5/gammaQ*T00-gammaP/gammaQ*T2000

    return result

def terTheta(phi1,phi2,phi3,phi4,N,p,q,r,s):
    alpha1=phi1.alpha[p]
    alpha2=phi2.alpha[q]
    alpha3=phi3.alpha[r]
    alpha4=phi4.alpha[s]

    gammaP=alpha1+alpha2
    gammaQ=alpha3+alpha4

    zetaP=alpha1*alpha2/gammaP
    zetaQ=alpha3*alpha4/gammaQ

    xi=gammaP*gammaQ/(gammaP+gammaQ)*1.0
    
    P=(alpha1*phi1.center+alpha2*phi2.center)/gammaP
    Q=(alpha3*phi3.center+alpha4*phi4.center)/gammaQ
    
    AB2=(LA.norm(phi1.center-phi2.center))**2
    CD2=(LA.norm(phi3.center-phi4.center))**2
    PQ2=(LA.norm(P-Q))**2

    #Kab=exp(-1*zetaP*AB2)
    #Kcd=exp(-1*zetaQ*CD2)
    A=phi1.center
    B=phi2.center
    C=phi3.center
    D=phi4.center
    ABX2=(A[0]-B[0])**2
    ABY2=(A[1]-B[1])**2
    ABZ2=(A[2]-B[2])**2

    CDX2=(C[0]-D[0])**2
    CDY2=(C[1]-D[1])**2
    CDZ2=(C[2]-D[2])**2


    Kab=exp(-1*zetaP*ABX2)*exp(-1*zetaP*ABY2)*exp(-1*zetaP*ABZ2)
    Kcd=exp(-1*zetaQ*CDX2)*exp(-1*zetaQ*CDY2)*exp(-1*zetaQ*CDZ2)

    t=xi*PQ2

    if N==0:
        if t<1e-14:
            F=1.0
        else:
            F=sqrt(pi)/2/sqrt(t)*erf(sqrt(t))

        total=2*(pi**(2.5))/gammaP/gammaQ/sqrt(gammaP+gammaQ)*\
                Kab*Kcd*F

    elif N==1:
        if t<1e-17:
            F=1/3.0
        elif t<0.18:
            F=coleTAYLOR(t,6,N)
        elif t<19.35:
            F=mac(t,N)
        else:
            #F=sqrt(pi)*.25*sqrt(pi)/t**1.5
            F=.25*sqrt(pi)/t**1.5
         


        total=2*pi**(2.5)/gammaP/gammaQ/sqrt(gammaP+gammaQ)*\
                Kab*Kcd*F


    else:
        if t<1e-17:
            F=1/float(2*N+1)
        elif t<0.18:
            F=coleTAYLOR(t,6,N)
        elif t<19.35:
            F=mac(t,N)
        else:
            F=0.375*sqrt(pi)/t**2.5

        total=2*pi**(2.5)/gammaP/gammaQ/sqrt(gammaP+gammaQ)*\
                    Kab*Kcd*F

    return total
