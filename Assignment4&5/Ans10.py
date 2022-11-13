import numpy as np 
import math
from math import sin,cos 
import matplotlib.pyplot as plt
from scipy.integrate import odeint
#import roboticstoolbox as rtb
#from spatialmath import SE3

#defining animation parameters 
fps = 600
dt = 1/fps
#define the geometric parameters of the robot

a1 = 4
a2 = 5
mlink1 = 1
mlink2 = 1
mlink3 = 1
radiuslink3 = 0.1
mp3 = mlink3
mt3 = mlink3
mR1 = 0
mR2 = 0
d = a1/2
I1 = (1/3)*(a1**2)*(mlink1)
IR1 = 0
I2 = (1/12)*(a2**2)*(mlink2)
IR2 = 0
I3 = I2 + mlink2*(d**2)
I4 = I1 + I3 + mlink2*(a1**2)
I5 = mlink2*a1*d
I6 = I3
I3_ = 0.5*mlink3*(radiuslink3)**2
I10 = mt3*(a1**2 + a2**2)
I8 = mp3*(a1**2 + a2**2) + I3_ + I10 
I9 = mp3*(a2**2) + I3_
I11 = 2*mp3*(a2**2) + 2*I3_ + I10
I12 = mp3*a1*a2
I13 = mt3*a1*a2
I14 = I4 + I8
I15 = I5 + I13
I16 = I6 + I9 + 0.25*I10 
I17 = I3 + I11
I18 = I5 + 2*(I12 + I13) 

# First step will be to solve the inverse kinematics for the SCARA robot.
# That is to calculate the desired joint variables from the given cartesian coordinates. 

def solveinvSCARA(p):
    x = p[0]
    y = p[1]
    z = p[2]
    d3 = -z
    D = ((x)**2 + (y)**2 - a1**2 - a2**2)/(2*a1*a2)
    q2 = math.acos(D)
    q1 = math.atan2((y),(x)) - math.atan2((a2*math.sin(q2)),(a1 + a2*math.cos(q2)))    
    #q1 = math.degrees(q1)
    #q2 = math.degrees(q2)
    return(q1,q2,d3)

# Second step is to create a PI controller for each joint
# The input to the controller will be the desired joint variable,the actual joint variable,Kp and Ki. 
# Here Kp is the proportional gain and Ki is the integral gain.  
# The output of the controller will be the torque that will be given to the dynamics of the system. 

def PIcontroller(qd,q,prev_e,dt,Kp,Ki):
    e = qd-q
    edt = prev_e + e*dt
    Tau_c = Kp*e + Ki*edt
    return (Tau_c,edt)

# Third step will be to incorporate dynamics.
# Here the tau for each joint will be mapped to corresponding joint variables. 
# Then to calculate joint variables from the provided tau (in this case the controller torque)
# the ODE will be solved. 
'''
def create_robot(a1,a2):
    robot = rtb.DHRobot(
        [
            rtb.RevoluteDH(a=a1,alpha=0,d=0),
            rtb.RevoluteDH(a=a2,alpha=math.pi,d=0),
            rtb.PrismaticDH(a=0,alpha=0,theta=0)
        ], name="scara3dof")
    return robot
   
class scara3dof(rtb.DHRobot):

    def __init__(self):
        super().__init__(
        [
            rtb.RevoluteDH(a=a1,alpha=0,d=0),
            rtb.RevoluteDH(a=a2,alpha=math.pi,d=0),
            rtb.PrismaticDH(a=0,alpha=0,theta=0)
        ], name="scara3dof")

scara = scara3dof()
'''

def Scara_dynamics_model(q,tau):
    q1 = q[0]
    q2 = q[1]
    q3 = q[2]
    q1dot = q[3]
    q2dot = q[4]
    q3dot = q[5]

    H = np.array([[(I14+2*I12*cos(q1) + 2*I15*cos(q2)),(0.5*(I17 + I18*cos(q2))),0],[(0.5*(I17 + I18*cos(q2))), (I16 + 0.5*I13*cos(q2)), 0],[0,0,mp3]])
    h = np.array([[-2*I15*q1dot*q2dot*sin(q2)-0.5*I18*(q2dot**2)*sin(q2)],[I15*(q1dot**2)*sin(q2)-0.25*I13*(q2dot**2)*sin(q2)],[0]])

    output_1 = np.array([[q1dot],[q2dot],[q3dot]])
    output_2 = np.matmul((np.linalg.inv(H)), (tau-h))
    output = np.vstack((output_1,output_2))
    return output

# Fourth step will be to solve the forward kinematics for the actual joint variabes

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

# Fifth step will be to define a function that takes the input as two cartesian points and returns the trajectory
# in terms of the joint variables.

def trajectory_planner(p_init,p_final):
    q_init = solveinvSCARA(p_init)
    q_final = solveinvSCARA(p_final)

    t0 = 0
    tf = 1

    q10 = q_init[0]
    q10_ = 0
    q1f = q_final[0]
    q1f_ = 0 

    q20 = q_init[1]
    q20_ = 0
    q2f = q_final[1]
    q2f_ = 0

    q30 = q_init[2]
    q30_ = 0
    q3f = q_final[2]
    q3f_ = 0

    T = [[1,t0,t0**2,t0**3],[0,1,2*t0,3*(t0**2)],[1,tf,tf**2,tf**3],[0,1,2*tf,3*(tf**2)]]
    T = np.array(T)
    Q1 = [q10,q10_,q1f,q1f_]
    Q1 = np.array(Q1)
    Q2 = [q20,q20_,q2f,q2f_]
    Q2 = np.array(Q2)
    Q3 = [q30,q30_,q3f,q3f_]
    Q3 = np.array(Q3)

    A = np.linalg.inv(T).dot(Q1)
    B = np.linalg.inv(T).dot(Q2)
    C = np.linalg.inv(T).dot(Q3)

    def fq1(t,A) :
        q = A[0] + t*A[1] + (t**2)*A[2] + (t**3)*A[3]
        q_dot = A[1] + (2*t)*A[2] + (3*(t**2))*A[3]
        #q_ddot = 2*A[2] + 6*t*A[3]
        return [q,q_dot]

    def fq2(t,B) :
        q = B[0] + t*B[1] + (t**2)*B[2] + (t**3)*B[3]
        q_dot = B[1] + (2*t)*B[2] + (3*(t**2))*B[3]
        #q_ddot = 2*B[2] + 6*t*B[3]
        return [q,q_dot]

    def fq3(t,C) :
        q = C[0] + t*C[1] + (t**2)*C[2] + (t**3)*C[3]
        q_dot = C[1] + (2*t)*C[2] + (3*(t**2))*C[3]
        #q_ddot = 2*C[2] + 6*t*C[3]
        return [q,q_dot]

    del_t = tf - t0
    N = del_t*(fps)
    t_ = t0
    q0 = [q10,q20,q30,q10_,q20_,q30_]

    run = True
    a = 1
    #speed = int((1/fps)*1000)
    trajectory = [q0]
    while run:
        if a <= N:
            #q0 = [fq1(t_,A)[0], fq2(t_,B)[0], fq1(t_,A)[1], fq2(t_,B)[1]]
            t_ += dt
            [q1d, q1d_] = [fq1(t_,A)[0], fq1(t_,A)[1]]
            [q2d, q2d_] = [fq2(t_,B)[0], fq2(t_,B)[1]]
            [q3d, q3d_] = [fq3(t_,B)[0], fq3(t_,B)[1]]
            #q_v = odeint(model,q0,[0,dt])
            #q1 = q_v[1][0]
            #q2 = q_v[1][1]
            trajectory.append([q1d,q2d,q3d,q1d_,q2d_,q3d_])
            a += 1 
        else:
            run = False
    
    return trajectory

# Sixth step is to plot the desired and actual trajectory on a 3D plot
def plot(traj):
    ax = plt.axes(projection='3d')

    # Data for a three-dimensional line
    xline = []
    yline = []
    zline = []
    for i in range(len(traj)):
        xline.append(traj[i][0])
        yline.append(traj[i][1])
        zline.append(traj[i][2])
    ax.plot3D(xline, yline, zline, 'red')
    #plt.show()

#plot([[1,2,2],[2,2,2],[3,2,2],[4,2,2],[5,2,2],[6,2,2]])

# 7th which is the last step is to integrate all the functions defined earlier and simulate the robot.

p_init = Scara_solveforward([0,0,0])
p_final = Scara_solveforward([30,60,6])
#print(p_init,p_final)
trajectory = trajectory_planner(p_init,p_final)

trajectory_c = []
for i in range(len(trajectory)):
    trajectory_c.append(Scara_solveforward([trajectory[i][0],trajectory[i][1],trajectory[i][2]]))
plot(trajectory_c)

q0 = np.array([solveinvSCARA(p_init)[0],solveinvSCARA(p_init)[1],solveinvSCARA(p_init)[2],0,0,0])
qf = np.array([solveinvSCARA(p_final)[0],solveinvSCARA(p_final)[1],solveinvSCARA(p_final)[2],0,0,0])
q = q0

prev_e_1 = 0
Kp1 = 1
Ki1 = 1

prev_e_2 = 0
Kp2 = 1
Ki2 = 1

prev_e_3 = 0
Kp3 = 1
Ki3 = 1

p_actual_traj = [p_init]

for i in range(len(trajectory)):
    qd = [trajectory[i][0],trajectory[i][1],trajectory[i][2]]
    tau_c_1,e_dt_1 = PIcontroller(qd[0],q[0],prev_e_1,dt,Kp1,Ki1)
    tau_c_2,e_dt_2 = PIcontroller(qd[1],q[1],prev_e_2,dt,Kp2,Ki2)
    tau_c_3,e_dt_3 = PIcontroller(qd[2],q[2],prev_e_3,dt,Kp3,Ki3)
    prev_e_1 = e_dt_1
    prev_e_2 = e_dt_3
    prev_e_2 = e_dt_3
    tau_c = np.array([[tau_c_1],[tau_c_2],[tau_c_3]])
    #Scara_dynamics_model(trajectory[i],tau_c)
    q_actual_dot = Scara_dynamics_model(trajectory[i],tau_c)
    q_actual = q + np.array([q_actual_dot[0][0]*dt,q_actual_dot[1][0]*dt,q_actual_dot[2][0]*dt,q_actual_dot[3][0]*dt,q_actual_dot[4][0]*dt,q_actual_dot[5][0]*dt])
    q = q_actual
    p_actual = Scara_solveforward([q_actual[0],q_actual[1],q_actual[2]])
    p_actual_traj.append(p_actual)
    traj = p_actual
    ax = plt.axes(projection='3d')

    # Data for a three-dimensional line
    xline = []
    yline = []
    zline = []
    for i in range(len(traj)):
        xline.append(traj[i][0])
        yline.append(traj[i][1])
        zline.append(traj[i][2])
    ax.plot3D(xline, yline, zline, 'red')
    plt.show()
#plot(p_actual_traj)
#plt.show()