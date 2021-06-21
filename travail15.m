%% by Bellila Ahmed Nassim, student at ENST
%On exchange at Université Laval
%Course UL GEL-7029 Prictive control under supervision of Prof. André Desbiens

%see Readme.txt to understand more the process


%%
clear all; clc;

%% pendulum parameters
K = 1.2;
m = 0.3;
L = 0.4;
g = 9.8;
% parameters – do not modify

%%
N = 150; % no. of simulation points, = 15 sec
T = 0.1; % sampling period 
Hc = 3;
Hp = 20;
lambda = 0.01;
R = 1; % Kalman - meas. noise var.
Q = 10*eye(3); % Kalman - process noise var.
P = 1*eye(3); % Kalman - P(0)
B = [0; T/(m*L^2)];
Ba = [0; T/(m*L^2); 0];

%% setpoint
r = [zeros(5,1); 450*ones(N-5,1)]; % in degrees
rrad = degtorad(r); % in radians
% disturbance
p = [zeros(50,1); 0.07*ones(N-50,1)];
% lower/upper bounds
lb = -1.5*ones(Hc,1);
ub = 1.5*ones(Hc,1);

%% initialisation
%initialize pendulum states
x = zeros(2,1);
% initialize estimated states
xe = zeros(3,1);
% initialize vector u
u = zeros(Hc,1);
% matrix to save the data (u and y)
simdata = zeros(N,2);

%% Simulation
for j = 1:N
 % Pendulum output
 y = x(1);
 
 % contol law
 % ukm1 = u(k-1):
 ukm1=u(1);
 % u = u(k):
 u = fmincon(@(u) evalcrit(Hp,Hc,T,lambda,xe,rrad(j),u,ukm1,K,m,L,g),...
 u,[],[],[],[],lb,ub) + p(j); 
 %u = u+.1;
 % Kalman filter
 F = [1 T 0; (-g*T/L)*cos(xe(1)) (K*T/m)+1 0; 0 0 1];
 H = [1 0 1];
 Kkalman = F*P*H'*inv(H*P*H'+R);
 %optimal states 
 xe = xe + Kkalman*(y-xe(1));
 xe = [(T*xe(2)+xe(1)); ((-g*T/L)*sin(xe(1))-(K*T/m)*xe(2)+xe(2)); xe(3)] + Ba*u(1);
 %update P
 P = F*(P-Kkalman*H*P)*F'+Q;
 % Pendulum states
 x = [T*x(2)+x(1); (-g*T/L)*sin(x(1))-(K*T/m)*x(2)+x(2)] + B*u(1);
 % Record data
 simdata(j,:) = [u(1) radtodeg(y)];
end


%% figures
figure(1);
t=(0:N-1)*T;
subplot(211);
stairs(t,simdata(:,1))
title('Manipulated variable')
subplot(212);
plot(t,simdata(:,2),'k',t,r,'r')
legend('output signal y','set-point r')
title('closed loop response')
xlabel('Time [s]');