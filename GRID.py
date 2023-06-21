from test import createMatrix
from test import clcTran

def gridFunction(NX, NY, NZ, DX, DY, DZ, LA, prosity, permability, DiActive):
    
 ### LT = Total Blocks  
    LT = NX*NY*NZ

 ### This Matix specifies active blocks and non-active blocks
 ### We can make ACTNUM functionality to better function our app
    ACTNUM = createMatrix(NX, NY, 0)
    if DiActive == True:
        for I in range(5,10): 
            ACTNUM[7][I] = -1
            ACTNUM[I][7] = -1

 ### Transformer matixes
    IT_IJ = createMatrix(NX, NY, 0)
    IJ_IT = createMatrix(LT, 2, 0)
    IA_IJ = createMatrix(NX, NY, 0)
    IJ_IA = createMatrix(LA, 2, 0)

    IT = 0
    IA = 0
    for IY in range(NY):
        
        for IX in range(NX-1,-1, -1):
            
            IT = IT+1
            IT_IJ[IX][IY] = IT
            IJ_IT[IT-1][0] = IX
            IJ_IT[IT-1][1] = IY
            

            if ACTNUM[IX][IY] == 0:
                IA = IA +1
                IA_IJ[IX][IY] = IA
                IJ_IA[IA-1][0] = IX
                IJ_IA[IA-1][1] = IY

    ### Rock  properties(prosity{fraction} permability{md}) nad pore volume
    PROS = createMatrix(LA, 1,  prosity)
    PERM = createMatrix(NX, NY, permability) 
    Prosity = createMatrix(LA, 1, 0.1)
    PV = createMatrix(LA, 1, 0) 
    for i in range(LA):
        PV[i][0]  = DX * DY * DZ * Prosity[i][0]


    ### Neighbour list (Clockwise approch for i & j )
    ### Transmissbilty matrix
    IY = 0
    IX = 0
    I = 0
    NL = createMatrix(LA, 4, 0)
    TL = createMatrix(LA, 4, 0)

    for IY in range(NY):
        for IX in range(NX-1, -1, -1):
            I = IA_IJ[IX][IY] 

            ### Up side of block(IX - 1)
            if  IX-1 == -1 or ACTNUM[IX-1][IY] == -1 :
                TL[I-1][0] = 0 
                NL[I-1][0] = -1
            else:
                NL[I-1][0] = IA_IJ[IX-1][IY]
                TL[I-1][0] = clcTran(DY, DZ, DY, PERM[IX][IY], PERM[IX-1][IY])
                
            ### Right side of block(IY + 1)
            if IY + 1 == 15 or ACTNUM[IX][IY+1] == -1  : 
                NL[I-1][1] = -1
                TL[I-1][1] = 0 
            else:
                NL[I-1][1] = IA_IJ[IX][IY+1]
                TL[I-1][1] = clcTran(DX, DZ, DY, PERM[IX][IY], PERM[IX][IY+1]) 

            ### Down side of block(IX - 1)
            if IX +1 == 15 or ACTNUM[IX+1][IY] == -1   : 
                NL[I-1][2] = -1
                TL[I-1][2] = 0
            else:
                NL[I-1][2] = IA_IJ[IX+1][IY]
                TL[I-1][2] = clcTran(DY, DZ, DY, PERM[IX][IY], PERM[IX+1][IY])

            ### Left side of block(IY - 1)
            if  IY-1 == -1 or ACTNUM[IX][IY-1] == -1 :  
                NL[I-1][3] = -1
                TL[I-1][3] = 0
            else:
                NL[I-1][3] = IA_IJ[IX][IY-1]
                TL[I-1][3] = clcTran(DX, DZ, DY, PERM[IX][IY], PERM[IX][IY-1])  

    return IT_IJ, IJ_IT, IA_IJ, IJ_IA, PV, NL, TL, LA, ACTNUM, PROS, PERM         