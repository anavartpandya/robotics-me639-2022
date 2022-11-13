import math
import numpy as np 

def solveinvSCARA(p,a1,a2):
    x = p[0]
    y = p[1]
    z = p[2]
    d3 = -z
    D = ((x)**2 + (y)**2 - a1**2 - a2**2)/(2*a1*a2)
    q2 = math.acos(D)
    q1 = math.atan2((y),(x)) - math.atan2((a2*math.sin(q2)),(a1 + a2*math.cos(q2)))    
    q1 = math.degrees(q1)
    q2 = math.degrees(q2)
    return(q1,q2,d3)

p = [0.8660254,1.5,-5] 
a1 = 1
a2 = 1 

ans = solveinvSCARA(p,a1,a2)
print(ans)