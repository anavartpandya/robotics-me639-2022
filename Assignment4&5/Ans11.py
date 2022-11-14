import numpy as np 
import math
from math import sin,cos 
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

######## SOME POINTS TO BE NOTED ########
# On running this code, there will be one matplotlib output window which will run a 3D animation.
# You will be able to see two lines, one green (which will be fixed) and the other will be red 
# which will be progressing further with time.
# The green line represents the desired trajectory and the red line represents the actual trajectory.
# The PI gains are tunned on the basis of current geometric parameters of the robot. So incase, you change the 
# geometric parameters, please also change the PI gains appropriately to achieve the desired trajectory. 

#defining animation parameters 

fps = 100
dt = 1/fps

#define the geometric parameters of the robot

a1 = 0
a2 = 4
a3 = 5

l1 = a1
lc1 = a1/2
l2 = a2
lc2 = a2/2
l3 = a3
lc3 = a3/2
m1 = 0
m2 = 1
m3 = 1

I1 = 0
I2 = (1/12)*m2*(l2**2)
I3 = (1/12)*m3*(l3**2)
I2_ = (1/24)*m2*(l2**2)
I3_ = (1/24)*m3*(l3**2)

# First step will be to solve the inverse kinematics for the PUMA robot.
# That is to calculate the desired joint variables from the given cartesian coordinates. 

def solveinvPUMA(p):
    x = p[0]
    y = p[1]
    z = p[2]
    q1 = math.atan2(y,x)
    w = (x**2 + y**2)**0.5    
    D = ((w)**2 + (z)**2 - a2**2 - a3**2)/(2*a2*a3)
    q3 = math.acos(D)
    q2 = math.atan2((z),(w)) - math.atan2((a3*math.sin(q3)),(a2 + a3*math.cos(q3)))    
    return(q1,q2,q3) 

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
# the ODE will be solved further for every time instance. 

def PUMA_dynamics_model(q,tau):
    q1 = q[0]
    q2 = q[1]
    q3 = q[2]
    q1dot = q[3]
    q2dot = q[4]
    q3dot = q[5]

    d11 = 0.5*(I2_ + I3_ + m2*(lc2*cos(q2))**2 + m3*(lc3*cos(q2+q3) + lc2*cos(q2))**2)
    d12 = 0
    d13 = 0
    
    d21 = 0
    d22 = m2*lc2**2 + m3*(l2**2 + lc3**2 + 2*l2*lc3*cos(q3)) + I2 + I3
    d23 = m3*(lc3**2 + l2*lc3*cos(q3)) + I3 

    d31 = 0
    d32 = m3*(lc3**2 + l2*lc3*cos(q3)) + I3
    d33 = m3*lc3**2 + I3

    Dmat = np.array([[d11,d12,d13],[d21,d22,d23],[d31,d32,d33]])

    h_ = -m3*l2*lc3*sin(q3)

    c111 = 0
    c121 = 0
    c131 = 0
    c211 = 0
    c221 = 0
    c231 = 0
    c311 = 0
    c321 = 0
    c331 = 0
    c1net = c111 + c121 + c131 + c211 + c221 + c231 + c311 + c321 + c331

    c112 = (-m2*0.5*(lc2*cos(q2)*sin(q2)) - m3*0.5*(lc3*cos(q2+q3) + l2*cos(q2))*(lc3*sin(q2+q3)+l2*sin(q2)))*q1dot*q1dot
    c122 = 0
    c132 = 0
    c212 = 0
    c222 = 0
    c232 = h_*q2dot*q3dot
    c312 = 0
    c322 = h_*q2dot*q3dot
    c332 = h_*q3dot*q3dot
    c2net = c112 + c122 + c132 + c212 + c222 + c232 + c312 + c322 + c332

    c113 = (-m3*0.5*(lc3*cos(q2+q3)+l2*cos(q2))*(lc3*sin(q2+q3)))*q1dot*q1dot
    c123 = 0
    c133 = 0
    c213 = 0
    c223 = -h_*q2dot*q2dot
    c233 = 0
    c313 = 0
    c323 = 0
    c333 = 0
    c3net = c113 + c123 + c133 + c213 + c223 + c233 + c313 + c323 + c333

    Cmat = np.array([[c1net],[c2net],[c3net]])

    output_1 = np.array([[q1dot],[q2dot],[q3dot]])
    output_2 = np.matmul((np.linalg.inv(Dmat)), (tau-Cmat))
    output = np.vstack((output_1,output_2))
    return output

# Fourth step will be to solve the forward kinematics for the actual joint variabes

def PUMA_solveforward(q): #input angles are in radians

    # taking input for joint space (in this case 3 angles)
    q1 = q[0]
    q2 = q[1]
    q3 = q[2]

    w = a2*cos(q2) + a3*cos(q2+q3)
    x = w*cos(q1)
    y = w*sin(q1)
    z = a2*sin(q2) + a3*sin(q2+q3)
    return [x,y,z]

# Fifth step will be to define a function that takes the input as two cartesian points and returns the trajectory
# in terms of the joint variables.

def trajectory_planner(p_init,p_final):
    q_init = solveinvPUMA(p_init)
    q_final = solveinvPUMA(p_final)

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
        return [q,q_dot]

    def fq2(t,B) :
        q = B[0] + t*B[1] + (t**2)*B[2] + (t**3)*B[3]
        q_dot = B[1] + (2*t)*B[2] + (3*(t**2))*B[3]
        return [q,q_dot]

    def fq3(t,C) :
        q = C[0] + t*C[1] + (t**2)*C[2] + (t**3)*C[3]
        q_dot = C[1] + (2*t)*C[2] + (3*(t**2))*C[3]
        return [q,q_dot]

    del_t = tf - t0
    N = del_t*(fps)
    t_ = t0
    q0 = [q10,q20,q30,q10_,q20_,q30_]

    run = True
    a = 1
    trajectory = [q0]
    while run:
        if a <= N:
            t_ += dt
            [q1d, q1d_] = [fq1(t_,A)[0], fq1(t_,A)[1]]
            [q2d, q2d_] = [fq2(t_,B)[0], fq2(t_,B)[1]]
            [q3d, q3d_] = [fq3(t_,C)[0], fq3(t_,C)[1]]
            trajectory.append([q1d,q2d,q3d,q1d_,q2d_,q3d_])
            a += 1 
        else:
            run = False
    
    return trajectory

# sixth step which is the last step is to integrate all the functions defined earlier and simulate the robot.

p_init = PUMA_solveforward([math.radians(0),math.radians(0),math.radians(0)])
p_final = PUMA_solveforward([math.radians(30),math.radians(60),math.radians(20)])
trajectory = trajectory_planner(p_init,p_final)

trajectory_c = []
for i in range(len(trajectory)):
    trajectory_c.append(PUMA_solveforward([trajectory[i][0],trajectory[i][1],trajectory[i][2]]))

q0 = np.array([solveinvPUMA(p_init)[0],solveinvPUMA(p_init)[1],solveinvPUMA(p_init)[2],0,0,0])
qf = np.array([solveinvPUMA(p_final)[0],solveinvPUMA(p_final)[1],solveinvPUMA(p_final)[2],0,0,0])
q = q0

prev_e_1 = 0
Kp1 = 295/2.5
Ki1 = 10

prev_e_2 = 0
Kp2 = 1500/2.5
Ki2 = 10

prev_e_3 = 0
Kp3 = 1000/2.5
Ki3 = 10

p_actual_traj = [p_init]

for i in range(len(trajectory)):
    qd = [trajectory[i][0],trajectory[i][1],trajectory[i][2]]
    e = qd[0] - q[0]
    tau_c_1,e_dt_1 = PIcontroller(qd[0],q[0],prev_e_1,dt,Kp1,Ki1)
    tau_c_2,e_dt_2 = PIcontroller(qd[1],q[1],prev_e_2,dt,Kp2,Ki2)
    tau_c_3,e_dt_3 = PIcontroller(qd[2],q[2],prev_e_3,dt,Kp3,Ki3)
    prev_e_1 = e_dt_1
    prev_e_2 = e_dt_3
    prev_e_2 = e_dt_3
    tau_c = np.array([[tau_c_1],[tau_c_2],[tau_c_3]])
    q_actual_dot = PUMA_dynamics_model(q,tau_c)
    q_actual = q + np.array([q_actual_dot[0][0]*dt,q_actual_dot[1][0]*dt,q_actual_dot[2][0]*dt,q_actual_dot[3][0]*dt,q_actual_dot[4][0]*dt,q_actual_dot[5][0]*dt])
    q = q_actual
    p_actual = PUMA_solveforward([q_actual[0],q_actual[1],q_actual[2]])
    p_actual_traj.append(p_actual)
    traj = p_actual_traj

plt.rcParams["figure.figsize"] = [7.00, 3.50]
plt.rcParams["figure.autolayout"] = True

def animate(num, data, line):
   line.set_alpha(0.7)
   line.set_data(data[0:2, :num])
   line.set_3d_properties(data[2, :num])
   return line

x_d = []
y_d = []
z_d = []
for i in trajectory_c:
    x_d.append(i[0])
    y_d.append(i[1])
    z_d.append(i[2])

x_a = []
y_a = []
z_a = []
for j in traj:
    x_a.append(j[0])
    y_a.append(j[1])
    z_a.append(j[2])

data = np.array([x_a, y_a, z_a])
N = len(traj)
fig = plt.figure()
ax = Axes3D(fig, auto_add_to_figure = False)
fig.add_axes(ax)

line, = plt.plot(data[0], data[1], data[2], lw=2, c='red')
ax.plot3D(x_d, y_d, z_d, 'green')
line_ani = animation.FuncAnimation(fig, animate, frames=N, fargs=(data, line), interval=50, blit=False)

plt.show()