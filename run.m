% initial conditions
clear all
close all

% Used Only for scalar systems

global K v N G A B n
A = [0 1; -1 0];
B = [0; 1];
n = length(A);
N = 3;
G = [0 1 1; 1 0 1; 1 0 0];
K1 = [-.4, -1];
K2 = [-.5, -.2];
K3 = [-.2, -.15];

K = [K1;K2;K3;];
v = [.13,.35];

% Initiallizing Parameters and Intial Conditions
n = length(A); %Size of system matrix
N = length(G); %Number of agents
% Initial conditions for each agent 
x10 = [-1 0]; % Agent 1
eta10 = 1;
tau10 = 0.1;
x20 = [1 0]; % Agent 2
eta20 = 2;
x30 = [0.5 0]; % Agent 3
eta30 = -1; 

x0 = [x10 x20 x30];
eta0 = [eta10 eta20 eta30];
tau0 = [tau10];
X0 = [x0 eta0 tau0]';
% simulation horizon
TSPAN=[0 20];
JSPAN = [0 1000];

% rule for jumps
% rule = 1 -> priority for jumps
% rule = 2 -> priority for flows
rule = 1;

options = odeset('RelTol',1e-6,'MaxStep',.1);

% simulate hybrid system using hybridsolver.m
[t, y, j] = hybridsolver( @f,@g,@C,@D,X0,TSPAN,JSPAN,rule,options,1);

% Define parameters P and sigma

P1 = [65.6696379735466,16.0354313062089,3.00912838963828,-8.65691900241049,8.12238678269490,-1.17184025397907,-10.0717918894850,2.14954950086325,-1.82855306036719;16.0354313062089,52.3761613383913,5.03160100549750,-13.2293606235898,-2.59482179554225,-1.83926221092715,-3.87977330149895,-5.66764358242275,-3.18166707949439;3.00912838963828,5.03160100549750,11.0182284281750,-1.45883022966958,-0.446375832283898,-4.84177026335474,-0.601026400978414,-0.257085741701844,-6.16422382578681;-8.65691900241049,-13.2293606235898,-1.45883022966958,67.3137207156636,2.47938184578875,6.53479664354068,-26.5507051994486,-7.65729499372210,-5.06964757035326;8.12238678269490,-2.59482179554225,-0.446375832283898,2.47938184578875,41.4036181148435,1.00981284729759,2.58906878331329,-17.3540176983263,-0.562570689440869;-1.17184025397907,-1.83926221092715,-4.84177026335474,6.53479664354068,1.00981284729759,14.9752706856879,-2.35045314033021,-1.01692421878622,-10.1215498716694;-10.0717918894850,-3.87977330149895,-0.601026400978414,-26.5507051994486,2.58906878331329,-2.35045314033021,14.6960585383124,1.59699413172436,2.95391114371741;2.14954950086325,-5.66764358242275,-0.257085741701844,-7.65729499372210,-17.3540176983263,-1.01692421878622,1.59699413172436,10.1026164345951,1.27533631768320;-1.82855306036719,-3.18166707949439,-6.16422382578681,-5.06964757035326,-0.562570689440869,-10.1215498716694,2.95391114371741,1.27533631768320,16.2973289426737];
sigma = 0.1;

% Define barAf, and barAg

Af = [A B; 0 0 0];
Bf = [zeros(2) B; 0 0 0];
Ag1 = [eye(2) zeros(2,1); K1, 0];
Ag2 = [eye(2) zeros(2,1); K2, 0];
Ag3 = [eye(2) zeros(2,1); K3, 0];

barAg = blkdiag(Ag1,Ag2,Ag3);
barAf = kron(eye(3),Af) - kron(diag(sum(G'))^-1*G,Bf);

% Evaluate Lyapunov function over the generated solutions
V = NaN(1,length(t));

% Building Error Matrix
Lp = eye(3) - diag(sum(G'))^-1*G;
for i = 1:length(t)
    err(:,i) = kron(Lp,eye(2))*y(i,[1:6])';
end

% Rearranging state vector
chi = [err([1,2],:); y(:,7)'; err([3,4],:); y(:,8)'; err([5,6],:); y(:,9)'];
tau = y(:,10);

for i =1:length(t)
    tau = y(:,N*n+N+1);
    V(i) = exp(sigma*tau(i))*chi(:,i)'*expm(barAf'*tau(i))*P1*expm(barAf*tau(i))*chi(:,i);
end

% Generate plots of data

x11 = y(:,1);
x12 = y(:,2);
x21 = y(:,3);
x22 = y(:,4);
x31 = y(:,5);
x32 = y(:,6);
figure(1)
set(1,'Position',[212 888 560 209])
subplot(1,2,2)
plot(t, V, 'Color', [.47 .67 .19]);
axis([0, 10, 0 240])
grid on
% axis([0,8,0,300])
subplot(1,2,1)
plot(x11,x12,'-.',x21,x22,':',x31,x32)
hold on
plot(x11(1),x12(1),'*',...
    'MarkerEdgeColor',[0,.45,.74])
plot(x21(1),x22(1),'*',...
    'MarkerEdgeColor',[0.85,.33,.1])
plot(x31(1),x32(1),'*',...
    'MarkerEdgeColor',[0.93,.69,.13])
grid on
legend(' ', ' ',' ')
%     
