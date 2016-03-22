function xdot = f(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab M-file                Author: Torstein Ingebrigtsen Bø
%
% Project: Simulation of a hybrid system (bouncing ball)
%
% Description: Flow map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global N A B n
n = length(A);
% differential equations
x1 = x([1,2]);
x2 = x([3,4]);
x3 = x([5,6]);
eta1 = x(7);
eta2 = x(8);
eta3 = x(9);
x1dot = A*x1 + B*eta1;
x2dot = A*x2 + B*eta2;
x3dot = A*x3 + B*eta3;

xdot = [x1dot; x2dot; x3dot];

etadot = zeros(N,1);
taudot = -1;
xdot = [xdot;etadot;taudot];

end