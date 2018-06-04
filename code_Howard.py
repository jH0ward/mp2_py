#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python
# Coleman Howard
# 5/9/11
#
# USAGE:
#   code_Howard.py > h2o.out &
#
# 
#
#
# MP2 Code from PHYS630
#   adapted from 661 HF code
#
#

import numpy as np
from numpy import zeros,zeros_like,shape
from numpy import linalg as LA
from numpy import matrix
from Gaussian import Gaussian3D
from basis import basis_set
from overlap import *
from kinetics import akinetic
import sys
from potential import ssAttraction
from potential import spAttraction
from potential import ppAttraction
from ter import ssssRepulsion
from ter import ssspRepulsion
from ter import ssppRepulsion
from ter import spspRepulsion
from reputable import ppppRepulsion
from reputable import pppsRepulsion
import normalize




##########################################################
#############    Check Input                 #############

infile='./h2o.in'
f=open(infile,'r')

### Begin gather info from input file
bname=f.readline().strip() # Grab basis set name
if bname=="STO3G":
    pass
elif bname ==  "ucSTO3G":
    pass
else:
    print "ERROR: Basis set not supported"
    sys.exit()

natom=int(f.readline().strip()) # Grab natom
atsym=[]
xyz=np.zeros(natom*3) # Define 3d cartesian array
xyz.shape=(natom,3)
Z=[]  # list that holds atomic charges
nucs=[] # list to hold x,y,z + atomic charge
bf = []

for i in range(natom):
    line_arry=f.readline().strip().split()
    atsym.append(line_arry[0])

    # Check that atomic symbol is supported by this program
    if atsym[i] == 'H':
        Z.append(1.0)
    elif atsym[i]=='He':
        Z.append(2.0)
    elif atsym[i]=='O':
        Z.append(8.0)
    else:
        sys.exit()
    
    xyz[i]=line_arry[1:] # Grab x,y,z for this atom

    nucs.append([])
    nucs[i].append(xyz[i])
    nucs[i].append(Z[i])

    # Define basis set for each atom
    basis_set(bf,atsym[i],bname,xyz[i])

for i in bf:
    l=i.l
    m=i.m
    n=i.n
    normalize.new(i)
    for j in range(i.nprims):
        res = i.NC*i.cc[j]*(2/pi)**0.75*(2**(l+m+n))\
                *i.alpha[j]**((2*l+2*m+2*n+3)/4.0)\
                /(dfactorial(2*l-1)*dfactorial(2*m-1)*dfactorial(2*n-1))**0.5
        i.cc[j]=res
## Let's try to normalize





f.close()
### End gather info from input file


######  Initialize Arrays
nbasis=len(bf)

##################################
## Basis Set Specification


####################################
####

#
#   _____  _____ ___   ___   ___  ___  ___
#     |      |   |__   |__|  |__|  |   |_
#   __|__    |   |__   |  \  |  |  |   |__
#

#### Begin iterations
S=np.zeros(nbasis**2)     # Overlap Matrix
S.shape=(nbasis,nbasis)    #  (n x n)
Hcore=np.zeros_like(S)   # Core hamiltonian holds 1-electron terms
Hcore.shape=(nbasis,nbasis)
G=np.zeros_like(S)       # 2-electron terms of Fock matrix 
F=np.zeros_like(S)       # Fock matrix holds Hcore + G
E=np.zeros_like(S)       # Diagonal E matrix holds orbital energies 
F2=np.zeros_like(S)      # Fock matrix in orthogonalized basis   
C=np.zeros_like(S)       # Coefficient matrix
C2=np.zeros_like(S)      # Coefficient matrix in orthogonalized basis
jk=np.zeros(nbasis**4)    # 2-electron integral matrix 
jk.shape=(nbasis,nbasis,nbasis,nbasis)



################################
#                              #
#     QC Integrals             #
#                              #
################################

S=np.zeros(nbasis**2)
S.shape=(nbasis,nbasis)

for i in range(len(bf)):
    nprims1=bf[i].nprims
    for j in range(i,len(bf)):
        nprims2=bf[j].nprims

        result=0.0
        for k in range(nprims1):
            cc1=bf[i].cc[k]

            for p in range(nprims2):
                cc2=bf[j].cc[p]

                result=result+(cc1*cc2*spoverlap(bf[i],bf[j],k,p))

        S[i][j]=result
        S[j][i]=S[i][j]

print "Overlap Integrals"
print S
# Compute kinetic integrals
print "\n\nComputing T Integrals"

for i in range(len(bf)):
    bf1=bf[i]
    nprims1=bf1.nprims
    for j in range(i,len(bf)):
        bf2=bf[j]
        nprims2=bf2.nprims

        result=0.0
        for k in range(nprims1):
            cc1=bf1.cc[k]

            for p in range(nprims2):
                cc2=bf2.cc[p]
                alpha2=bf2.alpha[p]

                result+=cc1*cc2*akinetic(bf1,bf2,k,p)

        Hcore[i][j]=result
        Hcore[j][i]=Hcore[i][j] # Take advantage of symmetry

print "\n\nKinetic Integrals"
print Hcore


# Compute nuclear attraction integrals
#   Here, we will update Hcore with nucler attraction (V) 
att=np.zeros_like(Hcore)
att.shape=(nbasis,nbasis)
print "\n\nComputing nuclear attraction integrals"
for i in range(len(bf)): # Loop over all (2) basis functions
    bf1=bf[i]
    nprims1=bf1.nprims

    for j in range(i,len(bf)):  # Start at i>j because Hcore is symm.
        bf2=bf[j]
        nprims2=bf2.nprims

        result=0.0
        for k in range(nprims1):
            cc1=bf1.cc[k]

            for p in range(nprims2):
                cc2=bf2.cc[p]

                if bf1.L + bf2.L == 0:
                    result+=cc1*cc2*ssAttraction(bf1,bf2,nucs,k,p)
                elif bf1.L + bf2.L == 1:
                    result+=cc1*cc2*spAttraction(bf1,bf2,nucs,k,p)
                elif bf1.L + bf2.L ==2:
                    result+=cc1*cc2*ppAttraction(bf1,bf2,nucs,k,p)
                
        #Hcore[i][j]+=result
        #Hcore[j][i]=Hcore[i][j]
        att[i][j]+=result
        att[j][i]=att[i][j]
print "Nuclear Attraction Integrals\n"
print att

Hcore+=att
# Compute electron repulsion integrals

# Determine unique orbital pairs
orb_pairs=[(i,j) for i in bf for j in bf if i.label<=j.label]

# This loop will perform calculations on unique pairs of unique pairs
count=0 
for i in range(len(orb_pairs)):  # This outer loop == left side of 2-electron integral
    bf1=orb_pairs[i][0]
    orb1 =bf1.label
    nprims1=bf1.nprims
    L1=bf1.L

    bf2=orb_pairs[i][1]
    orb2 = orb_pairs[i][1].label
    nprims2=bf2.nprims
    L2=bf2.L

    for j in range(i,len(orb_pairs)): # for j>=i, this loop is right side of 2-e int

        bf3=orb_pairs[j][0]
        orb3 = orb_pairs[j][0].label
        nprims3=bf3.nprims
        L3=bf3.L

        bf4=orb_pairs[j][1]
        orb4 = orb_pairs[j][1].label
        nprims4=bf4.nprims
        L4=bf4.L
        val=0.0
        whichint=[]


        # Create new send_pairs tuple
        sendME1=(bf1,bf2)
        sendME2=(bf3,bf4)

        for p in range(nprims1):
            cc1=bf1.cc[p]

            for q in range(nprims2):
                cc2=bf2.cc[q]

                for r in range(nprims3):
                    cc3=bf3.cc[r]

                    for s in range(nprims4):
                        cc4=bf4.cc[s]

                        which_prims = (p,q,r,s)

                        if L1+L2+L3+L4==0:
                            val+=ssssRepulsion(sendME1,sendME2,which_prims)*\
                                cc1*cc2*cc3*cc4
                        elif L1+L2+L3+L4==1:
                            val+=ssspRepulsion(sendME1,sendME2,which_prims)*\
                                cc1*cc2*cc3*cc4
                        elif L1+L2+L3+L4==4:
                            val+=ppppRepulsion(sendME1,sendME2,which_prims)*\
                                cc1*cc2*cc3*cc4
                        elif L1+L2+L3+L4==3:
                            val+=pppsRepulsion(sendME1,sendME2,which_prims)*\
                                cc1*cc2*cc3*cc4
                        elif L1+L2==2 or L3+L4==2:
                            val+=ssppRepulsion(sendME1,sendME2,which_prims)*\
                                cc1*cc2*cc3*cc4
                        elif L2+L4==2 or L1+L3==2 or L1+L4==2 or L2+L3==2:
                            val+=spspRepulsion(sendME1,sendME2,which_prims)*\
                                cc1*cc2*cc3*cc4
                           

        jk[orb1][orb2][orb3][orb4]=val

        # Impose 8-fold symmetry of 2-electron integrals
        jk[orb2][orb1][orb3][orb4]=jk[orb1][orb2][orb3][orb4]
        jk[orb2][orb1][orb4][orb3]=jk[orb1][orb2][orb3][orb4]
        jk[orb1][orb2][orb4][orb3]=jk[orb1][orb2][orb3][orb4]
        jk[orb3][orb4][orb1][orb2]=jk[orb1][orb2][orb3][orb4]
        jk[orb4][orb3][orb1][orb2]=jk[orb1][orb2][orb3][orb4]
        jk[orb4][orb3][orb2][orb1]=jk[orb1][orb2][orb3][orb4]
        jk[orb3][orb4][orb2][orb1]=jk[orb1][orb2][orb3][orb4]

print "2-e integrals"
print "2-e integral printing suppressed"
#print jk
#for i in range(nbasis):
#    for j in range(nbasis):
#      for k in range(nbasis):
#          for l in range(nbasis):
#              print i,j,k,l," ",jk[i][j][k][l]
iter=0        # 
P=np.zeros(nbasis**2) #Density matrix
P.shape=(nbasis,nbasis) #P is nbasis x nbasis
while iter<99:
    print "Beginning Iteration number ",iter
    print "-----------------------------"
    print "\n\n"
    print "\nComputing G matrix"
    # Add 2-electron integral terms to the empty G matrix
    G=np.zeros_like(S)
    G.shape=(nbasis,nbasis)
    for i in range(len(bf)):
        for j in range(len(bf)):
            for k in range(len(bf)):
                for l in range(len(bf)):
                    G[i][j]=G[i][j]+P[k][l]*(jk[i][j][k][l]-0.5*jk[i][l][k][j]) # Density matrix (P) scales 2-e contributions 
    print G

    # Add 1-electron contributions to form Fock matrix
    print "\nComputing Fock Matrix"
    F=Hcore+G
    print F

    ## Compute Electronic energy
    ee=0.0
    for i in range(len(bf)):
        for j in range(len(bf)):
            ee=ee+1.0/2.0*P[j][i]*(Hcore[i][j]+F[i][j])

    print "Electronic energy = ",ee

    evals,evecs=LA.eig(S)
    temp= matrix(S)*matrix(evecs)
    littleS=matrix(np.transpose(evecs))*temp
    littleS=array(littleS)


    half_s=np.zeros_like(littleS)
    for i in range(len(littleS)):
        half_s[i][i]=littleS[i][i]**(-1./2)


    temp = matrix(half_s)*matrix(np.transpose(evecs))
    bigS_half = matrix(evecs)*temp
    X=array(bigS_half)


    
    # Next line computes F'
    temp=matrix(F)*matrix(X)
    F2=matrix(np.transpose(X))*temp
    F2=np.array(F2)

    # Find eigenvalues of F' to find orbital energies
    evals,FX=LA.eig(F2)
    FXT=np.transpose(FX)

    # Must order orbital energies


    E=np.matrix(np.diag(evals,0)) # Add evals to E matrix
    Eun=np.copy(E)

    C2=FX
    sortC=np.zeros_like(C2)
    ## Find orbital energy order
    sortlist=np.sort(evals)
    docc=[]
    for i in range(nbasis):
        ind=np.searchsorted(sortlist,evals[i]) #search for evals[0]
        sortC[:,ind]=np.copy(C2[:,i])


    print "Sorted Orbital Energies"
    E=np.matrix(np.diag(sortlist,0))
    print E
    print "\nCoefficients Sorted"
    print sortC
    Cun=np.copy(C2)
    C2=np.copy(sortC)
    C=np.matrix(X)*matrix(C2) # Coefficients in original basis
    print "Updating Density Matrix"

    ## Update the density matrix
    oldP=P.copy()
    C=np.array(C)
    P=np.zeros_like(S)
    P.shape=(nbasis,nbasis)
    ndocc=5
    for i in range(nbasis):
        for j in range(nbasis):
            for k in range(ndocc):
                P[i][j]=P[i][j]+2*C[i][k]*C[j][k]
    
    print "Density Matrix"
    print P

    ## Compute root mean squaure
    print "Testing Convergence"
    rms=0.0
    for i in range(nbasis):
        for j in range(nbasis):
            rms=rms+(P[i][j]-oldP[i][j])**2
    rms=rms/nbasis**2
    print "Found rms value of ",rms

    ## Compute total energy
    R=1.4
    # Compute nuclear repulsion
    natom=len(nucs)
    nucRep=0.0
    for i in range(natom):
        for j in range(i+1,natom):
            nucRep+=nucs[i][1]*nucs[j][1]/LA.norm(nucs[i][0]-nucs[j][0])

    print "Total electronic Energy: ", ee
    print "Total Energy:  ",nucRep+ee
    
    if rms<1e-17:
        print "Converged Density Matrix in %d steps\n" % (iter+1)
        print "----------------------------------------------" 
        print "\tSCF Energy:\t",nucRep+ee
        print "----------------------------------------------" 
        break
    else:
        if iter>93:
            print "Maximum iterations reached: No convergence"
            sys.exit()
        print "Density matrix not converged:"
        print "Fill it up again, I'm gonna do one more"     
    iter=iter+1 
#
#               MP2
#   |
#
#
#
# Transform electron repulsion integrals
newjk=np.zeros_like(jk)
newjk.shape=(nbasis,nbasis,nbasis,nbasis)
for p in xrange(nbasis):
    for q in xrange(nbasis):
        for r in xrange(nbasis):
            for s in xrange(nbasis):
                val=0.0
                for i in xrange(len(C)):
                    Cip=C[i][p]
                    for j in xrange(len(C)):
                        Cjq=C[j][q]
                        for k in xrange(len(C)):
                            Ckr=C[k][r]
                            for l in xrange(len(C)):
                                Cls=C[l][s]
                                #newjk[p][q][r][s]+=C[i][p]*C[j][q]* \
                                # jk[i][j][k][l]*C[k][r]*C[l][s]
                                val+=Cip*Cjq* \
                                 jk[i][j][k][l]*Ckr*Cls
                newjk[p][q][r][s]+=val



# Compute mp2 energy
E=np.array(E)
emp2=0.0
docc=5
for i in range(ndocc):
    for a in range(ndocc,nbasis):
        for j in range(ndocc):
            for b in range(ndocc,nbasis):
                emp2+=newjk[i][a][j][b]*(2*newjk[i][a][j][b]-newjk[i][b][j][a])\
                   / (E[i][i]+E[j][j]-E[a][a]-E[b][b])

print "MP2 correlation E = ",emp2
print "----------------------------------------------" 
print "\tMP2 Total Energy:\t",emp2+nucRep+ee
print "----------------------------------------------" 
