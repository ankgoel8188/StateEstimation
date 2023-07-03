clc
clear all
close all
rand('state',1)
randn('state',1)
addpath './Functions/'
set(0, 'DefaultLineLineWidth', 1);

% System Definition
n = 3;
A = .3*rand(n,n);
B = rand(n,1);
C = rand(1,n);
Q = 2;1e-9;
R = 1;

a   = 0.8 + 0.6i;

Num = poly([0.3 0.2]);
Den = poly([a a' -0.5]);
% Den = poly([a a' 0.5]);
% G = tf(Num, Den,1)
% pzmap(G)
[A,B,C,D] = tf2ss(Num,Den);
B   = round(10*rand(size(A,2), 1))/10;
C   = round(10*randn(1,size(A,1)))/10;

% System Dimensions
lx  = size(A,1);
lu  = size(B,2);
lz  = size(C,1);
ly  = lz;

Q   = 1e-8*eye(lx);
R   = 1e-6*eye(ly);

% Initialize sizes
N       = 50;
x       = zeros(lx,N);
y       = zeros(ly,N);
u       = zeros(1,N);

xhat_OF = 0*x;

z       = zeros(ly,N);

P_OF    = 10*eye(n);
K_OF    = zeros(n, N);

J_OF    = zeros(1, N);

err_OF    = zeros(1, N);

x(:,1)  = [1 1 1]'+sqrt(Q)*randn(lx,1);
y(:,1)  = C*x(:,1)+sqrt(R)*randn(ly,1);


% Define structure to pass to filter function
System.A = A;
System.B = B;
System.C = C;
System.Q = Q;
System.R = R;


%% Simulation

for ii = 1:N
    
    % System simulation
    u(:,ii)     = 0*randn;
    x(:,ii+1)   = A*x(:,ii) + B*u(:,ii) + sqrt(Q)*randn(lx,1);
    y(:,ii+1)   = C*x(:,ii+1) + sqrt(R)*randn(ly,1);
    
    % Kalman Filter. 'OF' is the Optimal Filter (= Kalman Filter)
    [K_OF(:,ii),P_OF,J_OF(:,ii)] = Compute_KPJ(System, P_OF, 'OF' );
    xbar_OF         = A*xhat_OF(:,ii) + B*u(:,ii);
    xhat_OF(:,ii+1) = xbar_OF + K_OF(:,ii) * ( y(:,ii+1) - C*xbar_OF);
    err_OF(ii)      = norm(xhat_OF(:,ii)-x(:,ii))^2;
    
end


%% Plots
stairs(err_OF', 'linewidth',2)
hold on; axis tight; grid on
stairs(J_OF, 'linewidth',2)
set(gca,'yscale', 'log')
xlabel('Step')
legend('Error','Trace(P)')



