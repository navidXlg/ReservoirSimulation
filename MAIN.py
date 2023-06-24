from GRID import gridFunction
from test import X2CLCS
from test import EQU
from test import ACSAVE
from test import solution
from test import plot
import numpy as np
import copy


### Navid Khaleghi
### Reservoir NX = 15 , NY = 15 , NZ = 1;
### De-active = [(8 8), (8 9) (8 10) (8 7) (8 6) (9 8) (10 8) (7 8) (6 8) ]; 
### Block Size = { DX = 30m , DY = 30m , DZ = 10m };
### Prosity = 0.1 , Permabilty = 100 md;
### Reservoire Pressure = 4000psia , Sw = 0.1
### Simulation Time = 2 years  = 730day
#### Input to gridFunction ( NX, NY, NZ, DX, DY, DZ, LA, prosity, permability)

### 3.28 conversion factor (Meter to ft)
DX = 30 * 3.28
DY = 30 * 3.28
DZ = 10 * 3.28

### Number of Blocks
NX = 15
NY = 15
NZ = 1

### Prosity and Permability initial value
prosity = 0.1
permability = 100

### Number of Total Block
TA = 225

### Initial condition for X2CLC function 
Swc = 0.1 
Sor = 0.3

### Stability & Segregation for blocks P(PSIA) , saturation(fraction)
P1 = P2 =  4000
S1 = 0.1
S2 = 0.7
segregationIA = 119 

### Time duration of simulation and time steps (days)
Tmax = 730
Tstep = 10

### This function creates our grid system (True has deactive blocks - False does not have deactive blocks)
[IT_IJ, IJ_IT, IA_IJ, IJ_IA, PV, NL, TL, LA, ACTNUM, PROS, PERM] = gridFunction(NX, NY, NZ, DX, DY, DZ, TA, prosity, permability, False)

### Calling function for initialization
segrgatin = True
X = solution(TA, P1, S1, P2, S2, segrgatin, segregationIA )

### Plot of intial condition of reservoir
plot(X)

### Calling ACSAVE for next time step
SAVE = ACSAVE(X, LA)



### Calling x0 
Xo = copy.deepcopy(X)

### intial cahnge in X
if segrgatin == False:
     for IX in range(TA*2):
          X[IX][0] = X[IX][0] + 0.01*X[IX][0]

### Start of Simulation    
for dt in range(10, Tmax, Tstep):
     
     ### Error max  = 05
     sigma  = 0    

     ### Newton Raphson method
     for i in range(10):

          # Ft,Jt = EQU(X, X2, dt ,LA, ACTNUM, IJ_IT, NL, DX, DY, DZ, PROS, TL, SAVE[1], SAVE[0]) 
          Ft, Jt = EQU(X,SAVE[0] , SAVE[1] ,LA ,NL ,TL , dt, DX, DY, DZ)
          X = X - np.dot(np.linalg.inv(Jt), Ft)  #### Zarb matrix how to do it this x will replace initial guess 

          ### Clc the error for this NR
          for I in range (2*LA):
            sigma = sigma +abs(Xo[I][0] - X[I][0]) / Xo[I][0]
          Error = 1/(LA*2) * sigma  

          ### updating Xo
          for IX in range (2*LA):
             Xo[i][0] = X[i][0]

          if Error<10**-3:
               break


     ### X after my TimeStep
     print(X)
     
     ### Time OF X tuat we get after NR     
     print(dt)
     
     ### Plot for some time scales
     # tShow = [10, 100, 200, 300, 400, 500, 600, 700, 730]
     # if dt in tShow:
     plot(X)

     ### Calling ACSAVE for next time step
     SAVE = ACSAVE(X,LA)
     
          



        


