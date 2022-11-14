import numpy as np
import math

def solveinvSpherical(p,a1,a2):
    x = p[0]
    y = p[1]
    z = p[2]
    r = (x**2 + y**2)**0.5
    s = z - a1
    q1 = math.atan2(y,x)
    q2 = math.atan2(s,r)
    q3 = (r**2 + s**2)**0.5 - a2
    return [q1,q2,q3]

p = [0.8660254,1.5,-5] 
a1 = 1
a2 = 1 

ans = solveinvSpherical(p,a1,a2)
print(ans)

# Note that as we have not written code for forward kinematics of Spherical robot,
# we are not required to verify the results. This was mentioned by sir. 