import numpy as np
import math 
from math import sin, cos 
from operator import matmul
import pygame

height = 700
width = 700
window = pygame.display.set_mode((height,width))
backGroundColor=pygame.Color("WHITE")
window.fill(backGroundColor)
pygame_length_factor = int(300/1.5)


colour1 = (0,255,0)
colour2 = (255,0,0)
colour3 = (0,0,0)
circle_ = (350,350)
circle_radius = 10
border_width = 2

def to_pygame(coords, height):
    return (coords[0], height - coords[1]) 

l1 = 1*pygame_length_factor
l2 = 1*pygame_length_factor
l3 = 1*pygame_length_factor
radius = 1.5*pygame_length_factor
x1,y1 = 350,350

fps = 1000
dt = 1/fps
speed = 1

R = radius
T = 5
omega = 2*(math.pi/T)
t = 0
xd_list = []
yd_list = []
xdotd_list = []
ydotd_list = []
for i in range(int(T/dt)):
    xd = R*cos(omega*t)
    yd = R*sin(omega*t)
    xdotd = -R*omega*sin(omega*t)
    ydotd = R*omega*cos(omega*t)
    xd_list.append(xd)
    yd_list.append(yd)
    xdotd_list.append(xdotd)
    ydotd_list.append(ydotd)
    t += dt


def solve(N, DH_param, joints = []):  # DH = [a, alpha , d, theta], input alpha and theta in degrees.

    def T(N,DH):
        mat_list = []
        for i in range(N):
            a = DH[i][0]
            alpha = math.radians(DH[i][1])
            d = DH[i][2]
            theta = math.radians(DH[i][3])

            A_matrix = [[cos(theta), -(sin(theta))*cos(alpha), sin(theta)*sin(alpha), a*cos(theta)],[sin(theta), (cos(theta))*cos(alpha), -(cos(theta))*sin(alpha), a*sin(theta)], [0,sin(alpha),cos(alpha),d], [0,0,0,1]]
            A_matrix = np.array(A_matrix)
            mat_list.append(A_matrix)

        T_mat = mat_list[0]

        for j in mat_list[1:]:
            T_mat = matmul(T_mat, j)
        
        z_N = np.array([[T_mat[0][2]],[T_mat[1][2]],[T_mat[2][2]]])
        d_N = np.array([[T_mat[0][3]],[T_mat[1][3]],[T_mat[2][3]]])
        return (z_N,d_N)

    if joints == []: 
        joints = ['R']*N
    
    Jacob = None
    for i in range(len(joints)):
        if joints[i] == 'R':
            if i == 0:
                z_cross_O = np.transpose(np.cross(np.array([0,0,1]), np.transpose(T(N,DH_param)[1])))
                z = np.array([[0],[0],[1]])
                Jacob = np.concatenate((z_cross_O,z), axis = 0)
            else:
                z_cross_O = np.transpose(np.cross(np.transpose(T(i,DH_param)[0]), np.transpose(T(N,DH_param)[1] - T(i,DH_param)[1])))
                z = T(i,DH_param)[0]
                Jacob_i = np.concatenate((z_cross_O,z), axis = 0)
                Jacob = np.concatenate((Jacob,Jacob_i),axis = 1)    

        if joints[i] == 'P':
            if i == 0:
                z = np.array([[0],[0],[1]])
                zero_mat = np.array([[0],[0],[0]])
                Jacob = np.concatenate((z,zero_mat), axis = 0)
            else:
                z = T(i,DH_param)[0]
                zero_mat = np.array([[0],[0],[0]])
                Jacob_i = np.concatenate((z,zero_mat), axis = 0)
                Jacob = np.concatenate((Jacob,Jacob_i),axis = 1)
    Jacob = Jacob.round(4)
    return Jacob

N = 3 
joints_ = ['R','R','R'] 

def find_Jvinv(q1,q2,q3):
    DH_p = np.array([[l1,0,0,math.degrees(q1)],[l2,0,0,math.degrees(q2)],[l3,0,0,math.degrees(q3)]])  
    J = solve(N,DH_p,joints=joints_)
    Jv = J[:3]
    Jv_inv = np.linalg.pinv(Jv)
    return Jv_inv

def forward_kinematics(coord0,q1,q2,q3):
    x1_ = coord0[0] + l1*cos(q1)
    y1_ = coord0[1] + l1*sin(q1) 
    x2_ = coord0[0] + l1*cos(q1) + l2*cos(q1+q2)
    y2_ = coord0[1] + l1*sin(q1) + l2*sin(q1+q2)
    x3_ = coord0[0] + l1*cos(q1) + l2*cos(q1+q2) + l3*cos(q1+q2+q3)
    y3_ = coord0[1] + l1*sin(q1) + l2*sin(q1+q2) + l3*sin(q1+q2+q3)
    return([(x1_,y1_),(x2_,y2_),(x3_,y3_)])

Kp = 5
q_actual = []
q_prev = [0,math.acos(0.25), -2*math.acos(0.25)]
q1 = q_prev[0]
q2 = q_prev[1]
q3 = q_prev[2]
i = 0
run = True
while run:
    if i == len(xd_list):
        i = 0
    coord1 = (x1,y1)
    coords = forward_kinematics(coord1,q1,q2,q3)

    xdot = np.array([[xdotd_list[i]],[ydotd_list[i]],[0]])
    Jv_inv = find_Jvinv(q_prev[0],q_prev[1],q_prev[2])
    
    qdot = matmul(Jv_inv,xdot)

    q1 = q_prev[0] + (qdot[0][0])*dt
    q2 = q_prev[1] + (qdot[1][0])*dt
    q3 = q_prev[2] + (qdot[2][0])*dt
    
    coord1 = to_pygame((x1,y1), height)

    coord2 = coords[0]
    coord3 = coords[1]
    coord4 = coords[2]

    coord2 = to_pygame(coord2, height) 
    coord3 = to_pygame(coord3, height)
    coord4 = to_pygame(coord4, height)
    pygame.time.delay(speed)
    window.fill(backGroundColor)
    pygame.draw.line(window, colour1, coord1, coord2, border_width)
    pygame.draw.line(window, colour2, coord2, coord3, border_width)
    pygame.draw.line(window, colour2, coord3, coord4, border_width)

    pygame.draw.circle(window, colour3, coord1, circle_radius, border_width)
    pygame.draw.circle(window, colour3, coord2, circle_radius, border_width)
    pygame.draw.circle(window, colour3, coord3, circle_radius, border_width)
    pygame.draw.circle(window, colour3, coord4, circle_radius, border_width)
    pygame.draw.circle(window, colour3, coord1, radius, 1)

    pygame.display.update()
    q_prev = [q1,q2,q3]
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            run =False
    i += 1
pygame.quit()