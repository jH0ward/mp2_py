from Gaussian import *
import os,sys

def basis_set(bf,symbol,bname,coords):

    basis=eval(bname+"(symbol)")
    for i in range(0,len(basis),2):
        if basis[i]=='S':
            lvals = ( (0,0,0), )
            mark=1
            for j in lvals:
                bf.append(Gaussian3D(symbol,coords,basis,j,len(bf),mark))

        elif basis[i]=='2S':
            lvals = ( (0,0,0), )
            mark=3
            for j in lvals:
                bf.append(Gaussian3D(symbol,coords,basis,j,len(bf),mark))

        elif basis[i]=='2P':
            lvals = ( (1,0,0), (0,1,0), (0,0,1) )
            mark=5
            for j in lvals:
                bf.append(Gaussian3D(symbol,coords,basis,j,len(bf),mark))
        else:
            sys.exit()


def STO3G(symbol):

    basis =  { 'H': ('S', [
                          (3.42525091,0.15432897),
                          (0.62391373,0.53532814),
                          (0.16885540,0.44463454)
                        ]
                    ),
              'He': ('S', [
                          (6.36242139,0.15432897),
                          (1.15892300,0.53532814),
                          (0.31364979,0.44463454)
                        ]
                    ),
# I think the number of contraction shells is len(basis)/2
      #  the even indices of basis are nprims for the shell
      #  the  odd indices of basis are lists defining each primitive

                        # 1s functions
               'O': ('S', [ (130.70932000,0.15432897),
                          ( 23.8088610,0.53532814),
                          (  6.44360830,0.44463454)
                        ],

                        # 2s functions
                     '2S',[ (5.0331513,-0.09996723),
                          (1.1695961, 0.39951283),
                          (0.3803890, 0.70011547)
                          ],

                        # 2p functions
                     '2P', [ (5.0331513,0.15591627),
                            (1.1695961,0.60768372),
                            (0.3803890,0.39195739)
                          ]

                    )
            }



    return basis[symbol]

def ucSTO3G(symbol):
    cc = {'H' : (1.0,1.0,1.0),
                
          'He': (1.0,1.0,1.0),
                
          'O' : (1.0,1.0,1.0, \
                 1.0,1.0,1.0,1.0,\
                 1.0,1.0,1.0,1.0, \
                 1.0,1.0,1.0,1.0 \
                )
                 }
    alpha = {'H': (0.16885540,0.62391373,3.42525091),
            'He': (0.31364979,1.15892300,6.36242139),
            'O' : (6.4436083,23.8088610,130.7093200, \
                    0.3803890,0.3803890,0.380390,0.380390, \
                    1.1695961,1.1695961,1.1695961,1.1695961, \
                    5.0331513, 5.0331513, 5.0331513, 5.0331513 
                  )

            }
    return cc[symbol],alpha[symbol]
