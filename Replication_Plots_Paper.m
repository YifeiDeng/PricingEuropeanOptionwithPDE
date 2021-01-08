%% Replication of plots in the paper

%% ===================clear environment====================================

clear
close all
clc   
clf

%% Figure 1 and additional plot

% Numerical outputs for T = 0, .2, .4, .6, .8 and 1 yr for explicit method

% Parameter value from the paper:
K = 100;
S_min = 70;
S_max = 130;
r = 0.12;
T = 1.0;
sigma = 0.10;
dt = 0.01;
ds = 1.50;
N = T/dt;

% Obtain the option value at each expiry
[S,V,tao] = Explicit_B_S(K,S_max,r,T,sigma,ds,dt,'CALL');

% plot
figure(1)
plot(S,V(:,1)','-g',...
     S,V(:,round(N * 1/5))','-y',...
     S,V(:,round(N * 2/5))','-r',...
     S,V(:,round(N * 3/5))','-b',...
     S,V(:,round(N * 4/5))','-r',...
     S,V(:,N+1)','-k')
xlabel('Stock price (S)')
ylabel('Option Price (V(S,t))')
xlim([S_min S_max])
ylim([0 45])
title('Numerical solution for European call at different time steps with K = 100,r = 0.12, \sigma = 0.1, T = 1')
legend('time remaining to expiration = 1 yr', ...
       'time remaining to expiration = 0.8', ...
       'time remaining to expiration = 0.6', ...
       'time remaining to expiration = 0.4', ...
       'time remaining to expiration = 0.2', ...
       'time remaining to expiration = 0')
   
% 3-D plot of the option value, V(S,t)
figure(2)
mesh(tao,S,V)
title('European Call Option value, V(S,t), within the Explicit Method')
ylabel('Stock price (S)')
xlabel('Time (t)')

%% Figure 2

% Numerical outputs for T = 0, .2, .4, .6, .8 and 1 yr for Exact solution

% Parameter value from the paper:
T = [1.0, 0.8, 0.6, 0.4, 0.2, 0.0];

% Obtain analytic solution
V_exact = Exact_B_S(K,S_max,r,T,sigma,ds,'CALL');

% plot
figure(3)
plot(S,V_exact(:,1),'-g',...
     S,V_exact(:,2),'-y',...
     S,V_exact(:,3),'-r',...
     S,V_exact(:,4),'-b',...
     S,V_exact(:,5),'-r',...
     S,V_exact(:,6),'-k')
xlabel('Stock price (S)')
ylabel('Option Price (V(S,t))')
xlim([S_min S_max])
ylim([0 45])
title('Analytic solution for European call at different time steps with K = 100, r = 0.12,\sigma = 0.1, T = 1.')
legend('time remaining to expiration = 1 yr', ...
       'time remaining to expiration = 0.8', ...
       'time remaining to expiration = 0.6', ...
       'time remaining to expiration = 0.4', ...
       'time remaining to expiration = 0.2', ...
       'time remaining to expiration = 0')
   
 
%% Figure 3

% Compare the explicit solution with the analytic solution at inital time
figure(4)
hold on
plot(S,V(:,1),'-or')
plot(S,V_exact(:,1),'--g')
hold off
xlabel('Stock price (S)')
ylabel('Option Price (V(S,t))')
xlim([S_min S_max])
ylim([0 45])
title('Analytic solution and Numerical solution at initial time (K = 100, r = 0.12, \sigma = 0.1, T = 1)')
legend('Numerical Solution','Exact solution')

%% Figure 4:find the relative error for the scheme, the computation of the relative error in L1-norm

% Parameter value from the paper:
K = 100;
S_max = 130;
r = 0.12;
sigma = 0.10;
dt = 0.0500;
ds = 3;
T = 0:dt:1.0;

% Solution for the Explicit
[S,V,tao] = Explicit_B_S(K,S_max,r,T(end),sigma,ds,dt,'CALL');
% Solution for the Analytic formula
V_exact = Exact_B_S(K,S_max,r,T,sigma,ds,'CALL');

% Reverse the order of V_exact
V_exact = fliplr(V_exact);

% Initilzed the error vector
L_1_error = zeros(length(T),1);

% Calculate the relative error over time t from 0 to 1
for i  = 1:length(T)
    L_1_error(i) = vecnorm(V_exact(:,i)-V(:,i),1) / vecnorm(V_exact(:,i),1);
end

%plot the time versus relative error curve
figure(5)
plot(T,L_1_error,'-p')
xlabel('time (t)')
ylabel('Relative Error in L-1 Norm')
legend('Relative error for \Delta t = 0.0500 and \Delta x = 3')

%% Figure 5: the error is decreasing with respect to the smaller discretization parameters ?t and ?S
%% shows convergence of the explicit algorithm

%case #1 ? t=0.02 and ? x=2

dt1 = 0.0200;
ds1 = 2;
T1 = 0:dt1:1.0;

% Solution for the Explicit
[S,V,tao] = Explicit_B_S(K,S_max,r,T1(end),sigma,ds1,dt1,'CALL');
% Solution for the Analytic formula
V_exact = Exact_B_S(K,S_max,r,T1,sigma,ds1,'CALL');

% Initilzed the error vector
L_1_error_1 = zeros(length(T1),1);

% Reverse the order of V_exact
V_exact = fliplr(V_exact);

% Calculate the relative error over time t from 0 to 1
for i  = 1:length(T1)
    L_1_error_1(i) = vecnorm(V_exact(:,i)-V(:,i),1) / vecnorm(V_exact(:,i),1);
end

% Since at the end 0/0 will yield nan, we change it to 0
L_1_error_1(end) = 0.0;

%case #2 ? t=0.0125 and ? x=1.50

dt2 = 0.0125;
ds2 = 1.50;
T2 = 0:dt2:1.0;

% Solution for the Explicit
[S,V,tao] = Explicit_B_S(K,S_max,r,T2(end),sigma,ds2,dt2,'CALL');
% Solution for the Analytic formula
V_exact = Exact_B_S(K,S_max,r,T2,sigma,ds2,'CALL');

% Initilzed the error vector
L_1_error_2 = zeros(length(T2),1);

% Reverse the order of V_exact
V_exact = fliplr(V_exact);

% Calculate the relative error over time t from 0 to 1
for i  = 1:length(T2)
    L_1_error_2(i) = vecnorm(V_exact(:,i)-V(:,i),1) / vecnorm(V_exact(:,i),1);
end

%case #3 ? t= 0.0083 and ? x=1.20
dt3 = 0.0083;
ds3 = 1.20;
T3 = 0:dt3:1.0;

% Solution for the Explicit
[S,V,tao] = Explicit_B_S(K,S_max,r,T3(end),sigma,ds3,dt3,'CALL');
% Solution for the Analytic formula
V_exact = Exact_B_S(K,S_max,r,T3,sigma,ds3,'CALL');

% Initilzed the error vector
L_1_error_3 = zeros(length(T3),1);

% Reverse the order of V_exact
V_exact = fliplr(V_exact);

% Calculate the relative error over time t from 0 to 1
for i  = 1:length(T3)
    L_1_error_3(i) = vecnorm(V_exact(:,i)-V(:,i),1) / vecnorm(V_exact(:,i),1);
end

figure(6)
hold on
p1 = plot(T1,L_1_error_1,'-b');
p2 = plot(T2,L_1_error_2, '-g');
p3 = plot(T3,L_1_error_3, '-r');
hold off
xlabel('time (t)')
ylabel('Relative Error in L-1 Norm')
h = [p1(1);p2;p3(1)];
legend(h, 'Relative error for \Delta t = 0.02 and \Delta x = 2',...
        'Relative error for \Delta t = 0.0125 and \Delta x = 1.50',...
        'Relative error for \Delta t = 0.0083 and \Delta x = 1.20')
    
    
