function xplus = g(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab M-file               Author: Torstein Ingebrigtsen Bø
%
% Project: Simulation of a hybrid system (bouncing ball)
%
% Description: Jump map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global K N G v n

xg = x(1:N*n);

Ji = sum(G');

for i = 1:N
    e = kron(G(:,i),eye(n))'*(kron(ones(N,1),eye(n))*x((i-1)*n+1:i*n) - xg);
    etaplus(i) = K(i,:)*e/Ji(i);
end

tauplus = v(2) - (v(2) - v(1))*rand(1);

xplus = [x(1:n*N);etaplus';tauplus];

end