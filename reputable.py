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
        if marker>8:
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

def ppppRepulsion(pair1,pair2,myprims):
    tphi1=pair1[0]
    tphi2=pair1[1]
    tphi3=pair2[0]
    tphi4=pair2[1]

    tp=myprims[0]
    tq=myprims[1]
    tr=myprims[2]
    ts=myprims[3]

    xxyy='false'
    xxxy='false'
    # Reorder for convenience
    if tphi1.am==tphi3.am and tphi2.am==tphi4.am:
        phi1=tphi1
        p=tp
        phi2=tphi2
        q=tq
        phi3=tphi3
        r=tr
        phi4=tphi4
        s=ts
    elif tphi1.am==tphi4.am and tphi2.am==tphi3.am:
        phi1=tphi1
        p=tp
        phi2=tphi2
        q=tq
        phi3=tphi4
        r=ts
        pih4=tphi3
        s=tr
    elif tphi1.am==tphi2.am and tphi3.am==tphi4.am:
        xxyy='true'
        phi1=tphi1
        p=tp
        phi2=tphi2
        q=tq
        phi3=tphi3
        r=tr
        phi4=tphi4
        s=ts

    elif tphi1.am==tphi2.am and tphi1.am==tphi3.am and\
            tphi1.am!=tphi4.am:
        # Done
        xxxy='true'
        phi1=tphi1
        p=tp
        phi2=tphi2
        q=tq
        phi3=tphi3
        r=tr
        phi4=tphi4
        s=ts

    elif tphi1.am==tphi3.am and tphi1.am==tphi4.am and\
            tphi1.am!=tphi2.am:
        # Done
        xxxy='true'
        phi1=tphi3
        p=tr
        phi2=tphi4
        q=ts
        phi3=tphi1
        r=tp
        phi4=tphi2
        s=tq

    elif tphi1.am==tphi4.am and tphi1.am==tphi2.am and\
            tphi1.am!=tphi3.am:
        # Done
        xxxy='true'
        phi1=tphi1
        p=tp
        phi2=tphi2
        q=tq
        phi3=tphi4
        r=ts
        phi4=tphi3
        s=tr

    elif tphi2.am==tphi1.am and tphi2.am==tphi3.am and\
            tphi2.am!=tphi4.am:
        # Done
        xxxy='true'
        phi1=tphi1
        p=tp
        phi2=tphi2
        q=tq
        phi3=tphi3
        r=tr
        phi4=tphi3
        s=ts

    elif tphi2.am==tphi1.am and tphi2.am==tphi4.am and\
            tphi2.am!=tphi3.am:
        xxxy='true'
        # Done
        phi1=tphi1
        p=tp
        phi2=tphi2
        q=tq
        phi3=tphi4
        r=ts
        phi4=tphi3
        s=tr

    elif tphi2.am==tphi3.am and tphi2.am==tphi4.am and\
            tphi2.am!=tphi1.am:
        xxxy='true'
        # Done
        phi1=tphi3
        p=tr
        phi2=tphi4
        q=ts
        phi3=tphi2
        r=tq
        phi4=tphi1
        s=tp

    elif tphi3.am==tphi1.am and tphi3.am==tphi2.am and\
            tphi3.am!=tphi4.am:
        xxxy='true'
        # Done
        phi1=tphi1
        p=tp
        phi2=tphi2
        q=tq
        phi3=tphi3
        r=tr
        phi4=tphi4
        s=ts


    elif tphi3.am==tphi2.am and tphi3.am==tphi4.am and\
            tphi3.am!=tphi1.am:
        xxxy='true'
        # Done
        phi1=tphi3
        p=tr
        phi2=tphi4
        q=ts
        phi3=tphi2
        r=tq
        phi4=tphi1
        s=tp

    elif tphi3.am==tphi1.am and tphi3.am==tphi4.am and\
            tphi3.am!=tphi2.am:
        xxxy='true'
        # Done
        phi1=tphi3
        p=tr
        phi2=tphi4
        q=ts
        phi3=tphi1
        r=tp
        phi4=tphi2
        s=tq


    elif tphi4.am==tphi2.am and tphi4.am==tphi3.am and\
            tphi4.am!=tphi1.am:
        xxxy='true'
        # Done
        phi1=tphi3
        p=tr
        phi2=tphi4
        q=ts
        phi3=tphi2
        r=tq
        phi4=tphi1
        s=tp
    elif tphi4.am==tphi1.am and tphi4.am==tphi3.am and\
            tphi4.am!=tphi2.am:
        xxxy='true'
        # Done
        phi1=tphi3
        p=tr
        phi2=tphi4
        q=ts
        phi3=tphi1
        r=tp
        phi4=tphi2
        s=tq
    elif tphi4.am==tphi1.am and tphi4.am==tphi2.am and\
            tphi4.am!=tphi3.am:
        xxxy='true'
        # Done
        phi1=tphi1
        p=tp
        phi2=tphi2
        q=tq
        phi3=tphi4
        r=ts
        phi4=tphi3
        s=tr
    


    else:
        return 0.0
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

    # Determine if phi2 and phi4 are x,y,or z
    if phi1.l>0:
        mycomp1=0
    elif phi1.m>0:
        mycomp1=1
    elif phi1.n>0:
        mycomp1=2
    
    if xxyy=='false' and xxxy=='false':
        if phi2.l>0:
            comp2=0
        elif phi2.m>0:
            comp2=1
        elif phi2.n>0:
            comp2=2

    elif xxyy=='true':
        if phi3.l>0:
            comp2=0
        elif phi3.m>0:
            comp2=1
        elif phi3.n>0:
            comp2=2

    elif xxxy=='true':
        if phi4.l>0:
            comp2=0
        elif phi4.m>0:
            comp2=1
        elif phi4.n>0:
            comp2=2


    XPA=P[mycomp1]-A[mycomp1]
    XPB=P[mycomp1]-B[mycomp1]
    XPD=P[mycomp1]-D[mycomp1]
    XPQ=P[mycomp1]-Q[mycomp1]

    XCD=C[mycomp1]-D[mycomp1]

    XQP=Q[mycomp1]-P[mycomp1]
    XQD=Q[mycomp1]-D[mycomp1]
    XQC=Q[mycomp1]-C[mycomp1]
    XAB=A[mycomp1]-B[mycomp1]

    T00=terTheta(phi1,phi2,phi3,phi4,0,p,q,r,s)
    T01=terTheta(phi1,phi2,phi3,phi4,1,p,q,r,s)
    T02=terTheta(phi1,phi2,phi3,phi4,2,p,q,r,s)
    T03=terTheta(phi1,phi2,phi3,phi4,3,p,q,r,s)
    T04=terTheta(phi1,phi2,phi3,phi4,4,p,q,r,s)

    # Define T1
    T1=XPA*T00-xi/gammaP*XPQ*T01
    T1plus=XPA*T01-xi/gammaP*XPQ*T02
    T1plus2=XPA*T02-xi/gammaP*XPQ*T03
    T1plus3=XPA*T03-xi/gammaP*XPQ*T04

    # Define T2
    T2000=XPA*T1-xi/gammaP*XPQ*T1plus\
        +0.5/gammaP*(T00-xi/gammaP*T01)
    T2000plus=XPA*T1plus-xi/gammaP*XPQ*T1plus2\
        +0.5/gammaP*(T01-xi/gammaP*T02)

    T2000plus2=XPA*T1plus2-xi/gammaP*XPQ*T1plus3\
        +0.5/gammaP*(T02-xi/gammaP*T03)

    T2=(-1*alpha2*XAB-alpha4*XCD)/gammaQ*T1\
        +0.5/gammaQ*T00-gammaP/gammaQ*T2000

    T2plus=(-1*alpha2*XAB-alpha4*XCD)/gammaQ*T1plus\
        +0.5/gammaQ*T01-gammaP/gammaQ*T2000plus
    T2plus2=(-1*alpha2*XAB-alpha4*XCD)/gammaQ*T1plus2\
        +0.5/gammaQ*T02-gammaP/gammaQ*T2000plus2

    # Define T3
    T3000=XPA*T2000-xi/gammaP*XPQ*T2000plus\
        +1/gammaP*(T1-xi/gammaP*T1plus)

    T3000plus=XPA*T2000plus-xi/gammaP*XPQ*T2000plus2\
        +1/gammaP*(T1plus-xi/gammaP*T1plus2)


    T2010=(-1*alpha2*XAB-alpha4*XCD)/gammaQ*T2000\
          +1/gammaQ*T1-gammaP/gammaQ*T3000

    T2010plus=(-1*alpha2*XAB-alpha4*XCD)/gammaQ*T2000plus\
          +1/gammaQ*T1plus-gammaP/gammaQ*T3000plus

    T3=T2010+XAB*T2

    T3plus=T2010plus+XAB*T2plus

    # Define all T4 stuff
    T1100=T2000+XAB*T1
    T1100plus=T2000plus+XAB*T1plus
    T1100plus2=T2000plus2+XAB*T1plus2


    # Edit on talladega sunday
    #  i think i botched T0010
    T0010plus=XQC*T01-xi/gammaQ*XQP*T02
    T0010plus2=XQC*T02-xi/gammaQ*XQP*T03

    T0110plus=T2plus+XAB*T0010plus

    T1120=XQC*T3-xi/gammaQ*XQP*T3plus\
        +0.5/gammaQ*(T1100-xi/gammaQ*T1100plus)\
        +0.5/(gammaP+gammaQ)*(T0110plus+T2plus)

    T4=T1120+XCD*T3

    result=T4
    
    YQD=Q[comp2]-D[comp2]
    YPQ=P[comp2]-Q[comp2]
    YPB=P[comp2]-B[comp2]
    YQC=Q[comp2]-C[comp2]

    if xxxy=='true':
        result=YQD*T3+xi/gammaQ*YPQ*T3plus

    if phi1.am!=phi2.am:
        T1010=YQD*T2+xi/gammaQ*YPQ*T2plus
        T1010plus=YQD*T2plus+xi/gammaQ*YPQ*T2plus2
        result=YPB*T1010-xi/gammaP*YPQ*T1010plus+0.5/(gammaP+gammaQ)*\
                T2plus

    elif xxyy=='true':
        T11001=YQD*T1100+xi/gammaQ*YPQ*T1100plus
        T11001plus=YQD*T1100plus+xi/gammaQ*YPQ*T1100plus2

        result=YQC*T11001+xi/gammaQ*YPQ*T11001plus\
            +0.5/gammaQ*(T1100-xi/gammaQ*T1100plus)

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

    Kab=exp(-1*zetaP*AB2)
    Kcd=exp(-1*zetaQ*CD2)



    t=xi*PQ2

    if N==0:
        if t<1e-17:
            F=1.0
        else:
            F=sqrt(pi)/2/sqrt(t)*erf(sqrt(t))

        total=2*pi**(2.5)/gammaP/gammaQ/sqrt(gammaP+gammaQ)*\
                Kab*Kcd*F

    elif N==1:
        if t<1e-17:
            F=1/3.0
        elif t<0.18:
            F=coleTAYLOR(t,6,N)
        elif t<19.35:
            F=mac(t,N)
        else:
           # F=sqrt(pi)*.25*sqrt(pi)/t**1.5
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
            #F=sqrt(pi)*0.375*sqrt(pi)/t**2.5
            F=0.375*sqrt(pi)/t**2.5

        total=2*pi**(2.5)/gammaP/gammaQ/sqrt(gammaP+gammaQ)*\
                    Kab*Kcd*F

    return total

def pppsRepulsion(pair1,pair2,myprims):
# Switch the order

    tphi1=pair1[0]
    tphi2=pair1[1]
    tphi3=pair2[0]
    tphi4=pair2[1]

    #temp=phi4
    #phi4=phi1
    #phi1=temp

    tp=myprims[0]
    tq=myprims[1]
    tr=myprims[2]
    ts=myprims[3]
    if tphi1.L==0:
        phi4=tphi1
        s=tp
        phi1=tphi4
        p=ts
        phi3=tphi2
        r=tq
        phi2=tphi3
        q=tr

    elif tphi2.L==0:
        phi1=tphi3
        p=tr
        phi2=tphi4
        q=ts
        phi3=tphi1
        r=tp
        phi4=tphi2
        s=tq
    elif tphi3.L==0:
        phi1=tphi1
        p=tp
        phi2=tphi2
        q=tq
        phi3=tphi4
        r=ts
        phi4=tphi3
        s=tr

    else:
        phi1=tphi1
        p=tp
        phi2=tphi2
        q=tq
        phi3=tphi3
        r=tr
        phi4=tphi4
        s=ts

    
    xxys='false'
    xyxs='false'
    xyys='false'
    if phi1.am==phi2.am:
        if phi1.am!=phi3.am:
            xxys='true'
            if phi1.l==1:
                mycomp1=0
            elif phi1.m==1:
                mycomp1=1
            elif phi1.n==1:
                mycomp1=2
            if phi3.l==1:
                mycomp2=0
            elif phi3.m==1:
                mycomp2=1
            elif phi3.n==1:
                mycomp2=2

    elif phi1.am==phi3.am:
        if phi1.am!=phi2.am:
            #temp=phi2
            #tp=q
            #phi2=phi3
            #q=r
            #phi3=temp
            #r=tp
            xyxs='true'
            if phi1.l==1:
                mycomp1=0
            elif phi1.m==1:
                mycomp1=1
            elif phi1.n==1:
                mycomp1=2
            if phi2.l==1:
                mycomp2=0
            elif phi2.m==1:
                mycomp2=1
            elif phi2.n==1:
                mycomp2=2
    elif phi2.am==phi3.am:
        if phi2.am!=phi1.am:
            #temp=phi1
            #tp=p
            #phi1=phi3
            #p=r
            #phi3=temp
            #r=tp
            #xxys='true'
            xyys='true'

            if phi1.l==1:
                mycomp2=0
            elif phi1.m==1:
                mycomp2=1
            elif phi1.n==1:
                mycomp2=2
            if phi2.l==1:
                mycomp1=0
            elif phi2.m==1:
                mycomp1=1
            elif phi2.n==1:
                mycomp1=2

    if phi1.am!=phi2.am and phi1.am!=phi3.am and phi2.am!=phi3.am:
        return 0.0

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
    if xxys=='false' and xyys=='false' and xyxs=='false':
        if phi2.l>0:
            mycomp1=0
        elif phi2.m>0:
            mycomp1=1
        elif phi2.n>0:
            mycomp1=2

        if phi3.l>0:
            mycomp2=0
        elif phi3.m>0:
            mycomp2=1
        elif phi3.n>0:
            mycomp2=2

        if phi4.l>0:
            mycomp3=0
        elif phi4.m>0:
            mycomp3=1
        elif phi4.n>0:
            mycomp3=2

    XPA=P[mycomp1]-A[mycomp1]
    XPB=P[mycomp1]-B[mycomp1]
    XPD=P[mycomp1]-D[mycomp1]
    XPQ=P[mycomp1]-Q[mycomp1]

    XCD=C[mycomp1]-D[mycomp1]

    XQP=Q[mycomp1]-P[mycomp1]
    XQD=Q[mycomp1]-D[mycomp1]
    XQC=Q[mycomp1]-C[mycomp1]
    XAB=A[mycomp1]-B[mycomp1]

    T00=terTheta(phi1,phi2,phi3,phi4,0,p,q,r,s)
    T01=terTheta(phi1,phi2,phi3,phi4,1,p,q,r,s)
    T02=terTheta(phi1,phi2,phi3,phi4,2,p,q,r,s)
    T03=terTheta(phi1,phi2,phi3,phi4,3,p,q,r,s)
    T04=terTheta(phi1,phi2,phi3,phi4,4,p,q,r,s)

    # Define T1
    T1=XPA*T00-xi/gammaP*XPQ*T01
    T1plus=XPA*T01-xi/gammaP*XPQ*T02
    T1plus2=XPA*T02-xi/gammaP*XPQ*T03
    T1plus3=XPA*T03-xi/gammaP*XPQ*T04

    # Define T2
    T2000=XPA*T1-xi/gammaP*XPQ*T1plus\
        +0.5/gammaP*(T00-xi/gammaP*T01)
    T2000plus=XPA*T1plus-xi/gammaP*XPQ*T1plus2\
        +0.5/gammaP*(T01-xi/gammaP*T02)

    T2000plus2=XPA*T1plus2-xi/gammaP*XPQ*T1plus3\
        +0.5/gammaP*(T02-xi/gammaP*T03)

    T2=(-1*alpha2*XAB-alpha4*XCD)/gammaQ*T1\
        +0.5/gammaQ*T00-gammaP/gammaQ*T2000

    T2plus=(-1*alpha2*XAB-alpha4*XCD)/gammaQ*T1plus\
        +0.5/gammaQ*T01-gammaP/gammaQ*T2000plus

    # Define T3
    T3000=XPA*T2000-xi/gammaP*XPQ*T2000plus\
        +1/gammaP*(T1-xi/gammaP*T1plus)

    T3000plus=XPA*T2000plus-xi/gammaP*XPQ*T2000plus2\
        +1/gammaP*(T1plus-xi/gammaP*T1plus2)


    T2010=(-1*alpha2*XAB-alpha4*XCD)/gammaQ*T2000\
          +1/gammaQ*T1-gammaP/gammaQ*T3000

    T2010plus=(-1*alpha2*XAB-alpha4*XCD)/gammaQ*T2000plus\
          +1/gammaQ*T1plus-gammaP/gammaQ*T3000plus

    T3=T2010+XAB*T2

    if xxys=='true':
        YQC=Q[mycomp2]-C[mycomp2]
        YPB=P[mycomp2]-B[mycomp2]
        YPQ=P[mycomp2]-Q[mycomp2]
        T1100=T2000+XAB*T1
        T1100plus=T2000plus+XAB*T1plus
        T1000=T1
        T1000plus=T1plus
        T1000plus2=T1plus2
        T10001=YQC*T1000+xi/gammaQ*YPQ*T1000plus
        T10001plus=YQC*T1000plus+xi/gammaQ*YPQ*T1000plus2
        result=YQC*T1100+xi/gammaQ*YPQ*T1100plus
        return result

    if xyxs=='true':
        YPB=P[mycomp2]-B[mycomp2]
        YPQ=P[mycomp2]-Q[mycomp2]
        XQC=Q[mycomp1]-C[mycomp1]
        XPQ=P[mycomp1]-Q[mycomp1]
        T1000=T1
        T1000plus=T1plus
        result=YPB*T2-xi/gammaP*YPQ*T2plus#*T1010plus
        return result
    if xyys=='true':
        YPA=P[mycomp2]-A[mycomp2]
        YPQ=P[mycomp2]-Q[mycomp2]
        T1010=XQC*T1+xi/gammaQ*XPQ*T1plus\
            +0.5/(gammaP+gammaQ)*T01
        T1010plus=XQC*T1plus+xi/gammaQ*XPQ*T1plus2\
            +0.5/(gammaP+gammaQ)*T02
        XAB=A[mycomp1]-B[mycomp1]
        XQC=Q[mycomp1]-C[mycomp1]
        XQP=Q[mycomp1]-P[mycomp1]

        T00101=XQC*T01-xi/gammaQ*XQP*T02
        T00101plus=XQC*T02-xi/gammaQ*XQP*T03
        T0110=T1010+XAB*T00101
        T0110plus=T1010plus+XAB*T00101plus
        result=YPA*T0110-xi/gammaP*YPQ*T0110plus
        return result
    return T3
