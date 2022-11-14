import numpy as np
import math

# Here it is assumed that the Jacobian and the cartesian velocities are given.
# It is also assumed that the velocity input vector will be of nX1 dimension where n<=3. 
def solvejointspacevel(xdot,J):
    if len(J) == len(J[0]):
        qdot = np.linalg.solve(J,xdot)
        return qdot
    else:
        qdot = np.matmul(np.linalg.pinv(J),xdot) # for psuedo inverse cases
        return qdot

#example 

xdot = [1,3]
J = np.array([[1,2,3],[2,4,5]])
qdot = solvejointspacevel(xdot,J)
print(qdot)