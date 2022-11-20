import math 
import numpy as np 
#from math import cos, sin
import sympy as sp

t = sp.symbols('t')
#Defining Trajectory AB
t0 = 0
tf = 1

xAB = 0.4
zAB = -0.15

yAB0 = 0.06
yAB0_ = 0
yABf = 0.01
yABf_ = 0 

T = [[1,t0,t0**2,t0**3],[0,1,2*t0,3*(t0**2)],[1,tf,tf**2,tf**3],[0,1,2*tf,3*(tf**2)]]
T = np.array(T)
YAB = [yAB0,yAB0_,yABf,yABf_]
YAB = np.array(YAB)
A = np.linalg.inv(T).dot(YAB)
a0 = A[0]
a1 = A[1]
a2 = A[2]
a3 = A[3]

yAB = a0 + a1*t + a2*(t**2) + a3*(t**3)

# solving the inverse kinematics for joint variables calculation

a1 = 0.25
a2 = 0.25
d3max = 0.25
p = [xAB,yAB,zAB]

def solveinvSCARA(p):
    x = p[0]
    y = p[1]
    z = p[2]
    d3 = -z
    if abs(d3) > 0.25:
        return 'error'
    else:
        D = ((x)**2 + (y)**2 - a1**2 - a2**2)/(2*a1*a2)
        q2 = sp.acos(D)
        q1 = sp.atan2((y),(x)) - sp.atan2((a2*sp.sin(q2)),(a1 + a2*sp.cos(q2)))    
        return np.array([[q1],[q2],[d3]])

qAB = solveinvSCARA(p)
print('the caartesian coordinates as function of time are: ')
print('xAB = ')
print(xAB)
print()
print('yAB = ')
print(yAB)
print()
print('zAB = ')
print(zAB)
print()
print()
print('the joint variables as function of time are: ')
print('q1AB = ')
print(qAB[0][0])
print()
print('q2AB = ')
print(qAB[1][0])
print()
print('q3AB = ')
print(qAB[2][0])
print()
print('q1ABdot = ')
print(sp.diff(qAB[0][0]))
print()
print('q2ABdot = ')
print(sp.diff(qAB[1][0]))
print()
print('q3ABdot = ')
print(sp.diff(qAB[2][0]))
print()
print('q1ABddot = ')
print(sp.diff(sp.diff(qAB[0][0])))
print()
print('q2ABddot = ')
print(sp.diff(sp.diff(qAB[1][0])))
print()
print('q3ABddot = ')
print(sp.diff(sp.diff(qAB[2][0])))
print()

#Defining Trajectory BC
t0 = 0
tf = 1

yBC = 0.01
zBC = -0.15
xBC0 = 0.4
xBC0_ = 0
xBCf = 0.35
xBCf_ = 0 

T = [[1,t0,t0**2,t0**3],[0,1,2*t0,3*(t0**2)],[1,tf,tf**2,tf**3],[0,1,2*tf,3*(tf**2)]]
T = np.array(T)
XBC = [xBC0,xBC0_,xBCf,xBCf_]
XBC = np.array(XBC)
B = np.linalg.inv(T).dot(XBC)
b0 = B[0]
b1 = B[1]
b2 = B[2]
b3 = B[3]

xBC = b0 + b1*t + b2*(t**2) + b3*(t**3)

# solving the inverse kinematics for joint variables calculation
pBC = [xBC,yBC,zBC]

qBC = solveinvSCARA(pBC)
print('the caartesian coordinates as function of time are: ')
print('xBC = ')
print(xBC)
print()
print('yBC = ')
print(yBC)
print()
print('zBC = ')
print(zBC)
print()
print()
print('the joint variables as function of time are: ')
print('q1BC = ')
print(qBC[0][0])
print()
print('q2BC = ')
print(qBC[1][0])
print()
print('q3BC = ')
print(qBC[2][0])
print()
print('q1BCdot = ')
print(sp.diff(qBC[0][0]))
print()
print('q2BCdot = ')
print(sp.diff(qBC[1][0]))
print()
print('q3BCdot = ')
print(sp.diff(qBC[2][0]))
print()
print('q1BCddot = ')
print(sp.diff(sp.diff(qBC[0][0])))
print()
print('q2BCddot = ')
print(sp.diff(sp.diff(qBC[1][0])))
print()
print('q3BCddot = ')
print(sp.diff(sp.diff(qBC[2][0])))
print()

#Defining Trajectory CD
t0 = 0
tf = 1

xCD = 0.35
zCD = -0.15
yCD0 = 0.01
yCD0_ = 0
yCDf = 0.06
yCDf_ = 0 

T = [[1,t0,t0**2,t0**3],[0,1,2*t0,3*(t0**2)],[1,tf,tf**2,tf**3],[0,1,2*tf,3*(tf**2)]]
T = np.array(T)
YCD = [yCD0,yCD0_,yCDf,yCDf_]
YCD = np.array(YCD)
C = np.linalg.inv(T).dot(YCD)
c0 = C[0]
c1 = C[1]
c2 = C[2]
c3 = C[3]

yCD = c0 + c1*t + c2*(t**2) + c3*(t**3)

# solving the inverse kinematics for joint variables calculation
p = [xCD,yCD,zCD]

qCD = solveinvSCARA(p)
print('the caartesian coordinates as function of time are: ')
print('xCD = ')
print(xCD)
print()
print('yCD = ')
print(yCD)
print()
print('zCD = ')
print(zCD)
print()
print()
print('the joint variables as function of time are: ')
print('q1CD = ')
print(qCD[0][0])
print()
print('q2CD = ')
print(qCD[1][0])
print()
print('q3CD = ')
print(qCD[2][0])
print()
print('q1CDdot = ')
print(sp.diff(qCD[0][0]))
print()
print('q2CDdot = ')
print(sp.diff(qCD[1][0]))
print()
print('q3CDdot = ')
print(sp.diff(qCD[2][0]))
print()
print('q1CDddot = ')
print(sp.diff(sp.diff(qCD[0][0])))
print()
print('q2CDddot = ')
print(sp.diff(sp.diff(qCD[1][0])))
print()
print('q3CDddot = ')
print(sp.diff(sp.diff(qCD[2][0])))
print()
#Defining Trajectory DA
t0 = 0
tf = 1

yDA = 0.06
zDA = -0.15
xDA0 = 0.35
xDA0_ = 0
xDAf = 0.4
xDAf_ = 0 

T = [[1,t0,t0**2,t0**3],[0,1,2*t0,3*(t0**2)],[1,tf,tf**2,tf**3],[0,1,2*tf,3*(tf**2)]]
T = np.array(T)
XDA = [xDA0,xDA0_,xDAf,xDAf_]
XDA = np.array(XDA)
D = np.linalg.inv(T).dot(XDA)
d0 = D[0]
d1 = D[1]
d2 = D[2]
d3 = D[3]

xDA = d0 + d1*t + d2*(t**2) + d3*(t**3)

# solving the inverse kinematics for joint variables calculation
pDA = [xDA,yDA,zDA]

qDA = solveinvSCARA(pDA)
print('the caartesian coordinates as function of time are: ')
print('xDA = ')
print(xDA)
print()
print('yDA = ')
print(yDA)
print()
print('zDA = ')
print(zDA)
print()
print()
print('the joint variables as function of time are: ')
print('q1DA = ')
print(qDA[0][0])
print()
print('q2DA = ')
print(qDA[1][0])
print()
print('q3DA = ')
print(qDA[2][0])
print()
print('q1DAdot = ')
print(sp.diff(qDA[0][0]))
print()
print('q2DAdot = ')
print(sp.diff(qDA[1][0]))
print()
print('q3DAdot = ')
print(sp.diff(qDA[2][0]))
print()
print('q1DAddot = ')
print(sp.diff(sp.diff(qDA[0][0])))
print()
print('q2DAddot = ')
print(sp.diff(sp.diff(qDA[1][0])))
print()
print('q3DAddot = ')
print(sp.diff(sp.diff(qDA[2][0])))
print()