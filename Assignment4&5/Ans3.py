import numpy as np
import math

# Here it is assumed that the Jacobian and the cartesian velocities are given.
# It is also assumed that the velocity input vector will be of nX1 dimension where n<=3.
# It is also assumed that the input Jacobian will be a nXn square matrix (the part of the jacobian the 
# corresponds to the linear velocities term)
# To sum the above assumptions in one statement, I am assuming that the manipulator is not redundant. 
def solvejointspacevel(xdot,J):
    qdot = np.linalg.solve(J,xdot)
    return qdot

#example 

xdot = [1,3,4]
J = np.array([[1,2,3],[2,4,5],[3,-3,4]])
qdot = solvejointspacevel(xdot,J)
print(qdot)