function [qexp] = calc(qd)
qexp = sin(qd) + qd*qd;
