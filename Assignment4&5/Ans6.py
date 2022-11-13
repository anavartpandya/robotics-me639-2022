import math
import numpy as np

# Note that it is assumed that the input U (or R36) matrix is already provided.
# The function will output the euler angle corresponding to the spherical wrist. 
 
def solvespherical(U):
    u11 = U[0][0]
    u12 = U[0][1]
    u13 = U[0][2]
    u21 = U[1][0]
    u22 = U[1][1]
    u23 = U[1][2]
    u31 = U[2][0]
    u32 = U[2][1]
    u33 = U[2][2]
    if u13 != 0 and u23 != 0:
        theta1 = math.atan2(u33,(1-u33**2)**0.5)
        phi1 = math.atan2(u13,u23)
        psi1 = math.atan2(-u31,u32)

        theta2 = math.atan2(u33,-(1-u33**2)**0.5)
        phi2 = math.atan2(-u13,-u23)
        psi2 = math.atan2(u31,-u32) 
        return ([math.degrees(theta1),math.degrees(phi1),math.degrees(psi1)],[math.degrees(theta2),math.degrees(phi2),math.degrees(psi2)])

    elif u13 == 0 and u23 == 0: # infinitely many possible configurations

        theta1 = 0
        phi1 = 0 #by convention  
        psi1 = math.atan2(u11,u21) 

        theta2 = math.radians(180)
        phi2 = 0 #by convention  
        psi2 = -math.atan2(-u11,-u12)
        return ([math.degrees(theta1),math.degrees(phi1),math.degrees(psi1)],[math.degrees(theta2),math.degrees(phi2),math.degrees(psi2)])

#example
#input U (or R36) matrix
U = np.array([[1,2,3],[4,5,6],[7,8,1]])

euler_angles = solvespherical(U)
print('1st possible choice of euler angles :')
print(euler_angles[0])
print('2nd possible choice of euler angles :')
print(euler_angles[1])