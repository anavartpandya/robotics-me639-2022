function [qdot] = Robo2(tau,q)
m = 2;
c = 1;
b = 0.1;
q_ = q(1);
qdot_ = q(2);
qddot_ = (tau - c*qdot_ - b*q_ -sin(q_) -q_*q_)/m; 
qdot = [qdot_
        qddot_];