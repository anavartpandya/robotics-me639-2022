import numpy as np
import math
from math import sin,cos 

#Here I am using SCARA as the robot.
a1 = 0.25
a2 = 0.25
d3max = 0.25
#assuming the offset from the ground of the origin of SCARA to be 0.25 m, 
#I am setting the z coordinate as -0.25+0.1 = -0.15
pA = [0.45,0.075,-0.15]
pB = [0.45,-0.075,-0.15]
pC = [0.25,-0.075,-0.15]
pD = [0.25,0.075,-0.15]

def solveinvSCARA(p):
    x = p[0]
    y = p[1]
    z = p[2]
    d3 = -z
    if abs(d3) > 0.25:
        return 'error'
    else:
        D = ((x)**2 + (y)**2 - a1**2 - a2**2)/(2*a1*a2)
        q2 = math.acos(D)
        q1 = math.atan2((y),(x)) - math.atan2((a2*math.sin(q2)),(a1 + a2*math.cos(q2)))    
        return(math.degrees(q1),math.degrees(q2),d3)

def Scara_solveforward(q):

    def R(axis, angle): # Input axis as 'x', 'y' or 'z' strings for x, y or z axis respectively. Input angle in degrees
        angle = math.radians(angle)
        if axis == 'x':
            mat = np.array([[1,0,0],[0,math.cos(angle),-math.sin(angle)],[0,math.sin(angle),math.cos(angle)]])
        elif axis == 'y':
            mat = np.array([[math.cos(angle),0,math.sin(angle)],[0,1,0],[-math.sin(angle),0,math.cos(angle)]])
        elif axis == 'z':
            mat = np.array([[math.cos(angle),-math.sin(angle),0],[math.sin(angle),math.cos(angle),0],[0,0,1]])
        return mat

    #Defining H matrix calculation function
    def H_(R,d_):
        matH = np.zeros((4,4))
        for i in range(3):
            for j in range(3):
                matH[i][j] = R[i][j]
        for i in range(3):
            matH[i][3] = d_[i]
        matH[3][3] = 1
        return matH

    # taking input for joint space (in this case 2 angles and 1 displacement)
    q1 = q[0]
    q2 = q[1]
    d3 = q[2]

    # Calculating R matrices for SCARA 
    R01 = R('z',q1)
    R12 = R('z',q2) 
    R23 = R('z',0)

    # Calculating d vectors for SCARA
    d01 = np.array([0,0,0])
    d12 = np.array([a1,0,0])
    d23 = np.array([a2,0,0])

    # Calculating the Coordinate of the end effector with respect to the last coordinate frame 
    P3 = np.array([0,0,-d3])
    P3_ = np.array([P3[0],P3[1],P3[2],1])

    # Calculating H matrices for SCARA
    H01 = H_(R01,d01)
    H12 = H_(R12,d12)
    H23 = H_(R23,d23)

    # Calculating P0 Vector
    P0_ = np.matmul(np.matmul(H01,H12),np.matmul(H23,P3_))
    P0 = P0_[:3]

    return P0

# verifying point A
qA = solveinvSCARA(pA)
print('Joint variables for Point A are:')
print(qA)
P_A = Scara_solveforward(qA)
print('Ouput coordinate A is:')
print(P_A)

print()

# verifying point B
qB = solveinvSCARA(pB)
print('Joint variables for Point B are:')
print(qB)
P_B = Scara_solveforward(qB)
print('Ouput coordinate B is:')
print(P_B)

print()

# verifying point C
qC = solveinvSCARA(pC)
print('Joint variables for Point C are:')
print(qC)
P_C = Scara_solveforward(qC)
print('Ouput coordinate C is:')
print(P_C)

print()

# verifying point D
qD = solveinvSCARA(pD)
print('Joint variables for Point D are:')
print(qD)
P_D = Scara_solveforward(qD)
print('Ouput coordinate D is:')
print(P_D)
