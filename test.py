import numpy as np
import numpy as np
import matplotlib.pyplot as plt



### This function creates matrix based on row and colume and
### initial quntite
def createMatrix(row, colume, number):
    matrix = []
    for i in range(row):
        rowOFMatrix = []
        for j in range(colume):
            rowOFMatrix.append(number)
        matrix.append(rowOFMatrix)
    return matrix


### This function calcs Transmissibility based on half Transmissibility
def clcTran(Ax, Ay, Di, permi, permj):
    crossSection = Ax * Ay
    TransmissiI_J = ((permi * crossSection) / (Di/2)) * 0.001127
    TransmissiJ_I = ((permj * crossSection) / (Di/2)) * 0.001127
    Transmissibility = 1/((1/TransmissiJ_I)+(1/TransmissiI_J))

    return Transmissibility


### This function save special value
def ACSAVE(X, X2, LA):
    SAVEw = createMatrix(LA, 1, 0)
    SAVEo = createMatrix(LA, 1, 0)
    for i in range(LA):
        SAVEo[i][0] = (X2[i][18] *(1- X[2*i+1][0]))/X2[i][6]
        SAVEw[i][0] = (X2[i][18] *(X[2*i+1][0]))/X2[i][8]

    return SAVEo, SAVEw    


### This function clc the potential for block and its neighboure
def PotFunction(X, LA, X2):
    
    ### Gravity how to eork eith potential
    ### What is the pressuer of water
    g = 0
    Z = 0
    Zref = 15
    ### Z refrence for every block is 1000 meter
    ### depth of reservoir is 1000 meter
    ### Pot of oil and water
    POTo = createMatrix(LA, 1, 0)
    POTw = createMatrix(LA, 1, 0)

    for i in range(LA):
        POTo[i][0] = round(X[2*i][0], 6) - (X2[i][14] * Z)
        POTw[i][0] = round(X[2*i][0], 6) - X2[i][0] - (X2[i][15] * Z)  

    return POTo, POTw



                         
### This function gets primary varibels and returns secoundry varibels
def X2CLCS(X, Swc, Sor,  LT, PROS):

    ### Secondry varibels that we get from X1
    X2 = createMatrix(LT, 19, 0)
    Ros = 37.457 ###Ib/ft^3
    RS = 1000 ###SCF/bbl
    Rgs = 0.06248 ### Ib/ft^3
    Rowref = 62.3664 ### Ib/ft^3



    for i in range(LT):
        
        
        if (X[i*2+1][0]>1-Sor):
            X[i*2+1][0] = 1-Sor
            
        if(X[i*2+1][0]<Swc):
            X[i*2+1][0] = Swc
            
            
        Sr = 1/(1- Swc - Sor)
        Snw = (X[2*i+1][0] - Swc)/(1- Swc - Sor) 
        X2[i][0] = 0.34 * ((1 - Snw)**2) * 14.7  ### PC (PSIA)
        X2[i][1] = -2*(0.34)*(1-Snw)*Sr * 14.7    ### dPcSw (PSIA)  
        X2[i][2] = 0.4 * Snw**1.2   #### Krw 
        X2[i][3] = 0.48 * (Snw**0.2)*Sr #### dKrwSw
        X2[i][4] = 0.35 * (1 - Snw)**2.5  #### Kro
        X2[i][5] = - 0.35 * 2.5 * (1 - Snw)**1.5 * Sr ####dKroSw
        X2[i][6] = -6.0E-6 *X[2*i][0] + 1.275 ####Bo (bbl/STB)
        X2[i][7] = -6.0E-6     #####dBoP  (bbl/STB) d(PSIA)
        X2[i][8]= 1.03 * (1 + 3.0E-6 *(X[2*i][0] - 4000) + ((3.0E-6**2)*(X[2*i][0] - 4000) )/2) ####Bw (bbl/STB)
        X2[i][9] = 1.03 * (3E-6 + 9E-12 *(X[2*i][0] - 4000))    ##### dBwP (bbl/STB) d(PSIA)
        if (X[2*i][0]<2500):
            X2[i][10] = 0.5  ####MUo (cp)
            X2[i][11] = 0 ### dMUo (cp)
        elif (2500<X[2*i][0]<3000):
            X2[i][10] = 0.0001*X[2*i][0] + 0.25  
            X2[i][11] = 0.0001
        else:
            X2[i][10] = 0.6
            X2[i][11] = 0
        X2[i][12] = 0.4 #### MUw (cp)
        X2[i][13] = 0 ###dMUw   (cp)
        X2[i][14] = (Ros + ((RS * Rgs)/5.615)) / X2[i][6]   #### Roo (Ib/ft^3) errr
        X2[i][15] = Rowref / X2[i][8]  ###Row   (Ib/ft^3)
        X2[i][16] = X2[i][12] * X2[i][9] + X2[i][8] * X2[i][13]  ####dMUwBw = MUw * dBwP + Bw * dMUw
        X2[i][17] = X2[i][10] * X2[i][7] + X2[i][6] * X2[i][11]  ####dMUoBo = MUo * dBoP + Bo * dMUo
        X2[i][18] = PROS[i][0] * (1-0.000001*(X[2*i][0] - 4000) + (((0.000001)**2)*((X[2*i][0]-4000)**2))/2)  

    return X2



def solution(LA, P1, S1, P2, S2, segrgatin, segregationIA ):
    
    X = createMatrix(LA*2, 1, 0) ### (total blocks)
    if segrgatin  == True:
        ########## Segregation for blocks P(PSIA) , Saturation(fraction)
        for IX in range(LA):
              X[2*IX][0] = P1
              X[2*IX+1][0] = S1
             ### half of the blocks with 0.1 saturatuion
              if (IX>segregationIA):
                   X[2*IX][0] = P2
                   X[2*IX+1][0] = S2
    else:
       for IX in range(LA):
             X[2*IX][0] = P1
             X[2*IX+1][0] = S1               
    
    return X



### This function plot the distbution of Sw in reservoir

def plot(X):

    X2 = np.zeros((225, 1))

    for i in range(225):
        X2[i][0] =X[2*i+1][0] 
    # Create a 225x1 matrix of random values between 0 and 1
    
    # Reshape the matrix into a 15x15 grid
    grid = np.reshape(X2, (15, 15))
    
    
    # Create a color map with 256 levels, ranging from blue to red
    cmap = plt.cm.get_cmap('RdBu', 256)
    
    # Normalize the values to fit within the color map range
    norm = plt.Normalize(vmin=X2.min(), vmax=X2.max())
    
    # Create a 15x15 grid of squares with colors based on the property values
    fig, ax = plt.subplots()
    for i in range(grid.shape[0]):
        for j in range(grid.shape[1]):
            color = cmap(norm(grid[i, j]))
            rect = plt.Rectangle((j, i), 1, 1, facecolor=color, edgecolor='black')
            ax.add_patch(rect)
    
    # Set the x and y limits to show the entire grid
    ax.set_xlim([0, grid.shape[1]])
    ax.set_ylim([0, grid.shape[0]])
    
    # Add a color bar to show the color map legend
    cbar = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax)
    cbar.set_label('Water Saturation')
    
    # Add axis labels and title
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Sw')
    
    # Show the plot
    plt.show()




######### EQU function clc the jacobin and redisule base on X2
def EQU (X, X2, ACSAVE_o, ACSAVE_w, LA, NL, TL, dt, DX, DY, DZ):
        
        ### Clc bulk volume for every block
        Vb = DX * DY * DZ
    
        ### Jacobin and Residul Matrixes
        F = createMatrix (LA*2,1,0)
        J = createMatrix (LA*2,LA*2,0)

        ### Potential Function for specifing whether block is upstrem or downstream
        [POT_o , POT_w] = PotFunction(X, LA, X2)
        
        for IA in range(LA):
            
            ### Structure matrix for every block to complete main matrixes
            JB = createMatrix(2,10,0)
            FB = createMatrix(2,1,0)

            for n in range (4):

                #### neighbors is active or inactive
                if (NL[IA][n] > -1):

                    JA = NL[IA][n]
                    JA = JA - 1


                    ### OIl Neighbor  is Upstream
                    if (POT_o[JA][0] >= POT_o[IA][0]):
                        JB[0][0] = JB[0][0] + TL[IA][n] * (X2[JA][4] / (X2[JA][10] * X2[JA][6])) * (-1) 
                        JB[0][1] = JB[0][1] + 0 
                        JB[0][2*n + 2] = TL[IA][n] * (-X2[JA][17]) * (X2[JA][4] / (X2[JA][6] * X2[JA][10])**2) * (POT_o[JA][0] - POT_o[IA][0]) + TL[IA][n] * (X2[JA][4] / (X2[JA][10] * X2[JA][6])) * (1) 
                        JB[0][2*n + 3] = TL[IA][n] * (1 / (X2[JA][6] * X2[JA][10])) * X2[JA][5] * (POT_o[JA][0] - POT_o[IA][0])
                        FB[0][0] += TL[IA][n] * (X2[JA][4] / (X2[JA][10] * X2[JA][6])) * (POT_o[JA][0] - POT_o[IA][0])
                    
                    ### oil Central  is Upstream
                    else:
                        JB[0][0] = JB[0][0] + TL[IA][n] * (-X2[IA][17]) * (X2[IA][4] / (X2[IA][6] * X2[IA][10])**2) * (POT_o[JA][0] - POT_o[IA][0]) + TL[IA][n] * (X2[IA][4] / (X2[IA][6] * X2[IA][10])) * (-1) 
                        JB[0][1] = JB[0][1] + TL[IA][n] * (1 / (X2[IA][6] * X2[IA][10])) * X2[IA][5] * (POT_o[JA][0] - POT_o[IA][0]) 
                        JB[0][2*n + 2] = TL[IA][n] * (X2[IA][4] / (X2[IA][6] * X2[IA][10] ))
                        JB[0][2*n + 3] = 0
                        FB[0][0] += TL[IA][n] * (X2[IA][4] / (X2[IA][10] * X2[IA][6])) * (POT_o[JA][0] - POT_o[IA][0])

                    
                    ### Water Neighbor  is Upstream 
                    if (POT_w[JA] >= POT_w[IA]):
                        JB[1][0] = JB[1][0] + TL[IA][n] * (X2[JA][2] / (X2[JA][12] * X2[JA][8])) * (-1)   
                        JB[1][1] = JB[1][1] + TL[IA][n] * (X2[JA][2] / (X2[JA][12] * X2[JA][8])) * (X2[IA][1]) 
                        JB[1][2*n + 2] = TL[IA][n] * (-X2[JA][16]) * (X2[JA][2] / (X2[JA][12] * X2[JA][8])**2) * (POT_w[JA][0] - POT_w[IA][0]) + TL[IA][n] * (X2[JA][2] / (X2[JA][12] * X2[JA][8])) * (1)
                        JB[1][2*n + 3] = TL[IA][n] * (1/ (X2[JA][12] * X2[JA][8])) * X2[JA][3] * (POT_w[JA][0] - POT_w[IA][0]) + TL[IA][n] * (X2[JA][2] / (X2[JA][12] * X2[JA][8])) * (-X2[JA][1])
                        FB[1][0] += TL[IA][n] * (X2[JA][2] / (X2[JA][12] * X2[JA][8])) * (POT_w[JA][0] - POT_w[IA][0]) 

                    ### Water center is Upstream
                    else:
                        JB[1][0] = JB[1][0] + TL[IA][n] * (X2[IA][2] / (X2[IA][12] * X2[IA][8])**2) * (-X2[IA][16]) * (POT_w[JA][0] - POT_w[IA][0]) + TL[IA][n] * (X2[IA][2] / (X2[IA][12] * X2[IA][8])) * (-1) 
                        JB[1][1] = JB[1][1] + TL[IA][n] * (1/ (X2[IA][12] * X2[IA][8])) * X2[IA][3] * (POT_w[JA][0] - POT_w[IA][0]) + TL[IA][n] * (X2[IA][2] / (X2[IA][12] * X2[IA][8])) * (X2[IA][1]) 
                        JB[1][2*n + 2] = TL[IA][n] * (X2[IA][2] / (X2[IA][12] * X2[IA][8])) * 1
                        JB[1][2*n + 3] = TL[IA][n] * (X2[IA][2] / (X2[IA][12] * X2[IA][8])) * (-X2[JA][1])
                        FB[1][0] += TL[IA][n] * (X2[IA][2] / (X2[IA][12] * X2[IA][8])) * (POT_w[JA][0] - POT_w[IA][0])



            ### Adding Accumulation term for every Block
            JB[0][0] += - Vb/dt  / 5.615 * ((-X2[IA][7]) * (X2[IA][18]*(1-X[2*IA+1][0]) / (X2[IA][6] **2))) 
            JB[0][1] += - Vb/dt  / 5.615 * (X2[IA][18] / X2[IA][6] * (-1)) 
            JB[1][0] += - Vb/dt  / 5.615 * ((-X2[IA][9]) * (X2[IA][18] * X[2*IA+1][0] / X2[IA][8] **2)) 
            JB[1][1] += - Vb/dt  / 5.615 *((X2[IA][18] / X2[IA][8]) * 1) 
            FB[0][0] += -Vb/dt  / 5.615 * ((X2[IA][18] * (1-X[2*IA+1][0]) / X2[IA][6]) - ACSAVE_o[IA][0] ) 
            FB[1][0] += -Vb/dt  / 5.615 *((X2[IA][18] * X[2*IA+1][0] / X2[IA][8]) - ACSAVE_w[IA][0] )

            ### Adding structer matrix to main matrixes
            F[2*IA][0] = FB[0][0]
            F[2*IA+1][0] = FB[1][0]
            J[2*IA][2*IA] = JB[0][0]
            J[2*IA][2*IA+1] = JB[0][1]
            J[2*IA+1][2*IA] = JB[1][0]
            J[2*IA+1][2*IA+1] = JB[1][1]
            for n in range (4):
                if (NL[IA][n] != -1):
                    JA = NL[IA][n]

                    J[2*IA][2*JA-2] = JB[0][2*n+2]
                    J[2*IA][2*JA-1] = JB[0][2*n+3]
                    J[2*IA+1][2*JA-2] = JB[1][2*n+2]
                    J[2*IA+1][2*JA-1] = JB[1][2*n+3]
                    
            

        return F,J


### Defing object for every block
class Block:

    def __init__(self, Swc, Sor, P, Sw):

        self.Swc = Swc
        self.Sor = Sor
        self.Sw = Sw
        self.P = P
        self.dBoP = -6.0E-6
        self.Sr = 1/(1- Swc - Sor)
        self.Snw =  (Sw - Swc)/(1- Swc - Sor) 
        self.MUw = 0.4
        self.dMUw = 0
        self.Ros = 37.457 ###Ib/ft^3
        self.RS = 1000 ###SCF/bbl
        self.Rgs = 0.06248 ### Ib/ft^3
        self.Rowref = 62.3664 ### Ib/ft^3


    def PC(self):
        return 0.34 * ((1 - self.Snw)**2) * 14.7
    
    def dPcSw(self):
        return -2*(0.34)*(1- self.Snw )* self.Sr * 14.7
     
    def Krw (self):
        return 0.4 * (self.Snw**1.2)
    
    def dKrwSw(self):
        return 0.48 * (self.Snw**0.2)*self.Sr 
    
    def  Kro(self):
        return 0.35 * (1 - self.Snw)**2.5 
    
    def dKroSw(self):
        return - 0.35 * 2.5 * (1 - self.Snw)**1.5 * self.Sr
    
    def Bo(self):
        return -6.0E-6 *self.P + 1.275
    
    def Bw(self):
        return 1.03 * (1 + 3.0E-6 *(self.P - 4000) + ((3.0E-6**2)*(self.P - 4000) )/2) 
    
    def dBwP(self):
        return 1.03 * (3E-6 + 9E-12 *(self.P - 4000))
    
    def MUo(self):

        if (self.P<2500):
            MUo = 0.5  ####MUo (cp)
        elif (2500<self.P<3000):
            MUo = 0.0001*self.P + 0.25  
        else:
            MUo = 0.6
        return MUo 

    def dMUo (self):

        if (self.P<2500):
            dMUo = 0 ### dMUo (cp)
        elif (2500<self.P<3000): 
            dMUo = 0.0001
        else:
            dMUo = 0
        return dMUo   

    def Roo(self):
        return (self.Ros + ((self.RS * self.Rgs)/5.615)) / self.Bo()

    def Row(self):
        return self.Rowref / self.Bw()
    
    def dMUwBw (self):
        return self.MUw() * self.dBwP() + self.Bw() * self.dMUw
    
    def dMUoBo(self):
        return self.MUo() * self.dBoP + self.Bo() * self.dMUo()

    def PROS(self):
        return 0.1 * (1-0.000001*(self.P - 4000) + (((0.000001)**2)*((self.P-4000)**2))/2)  
    