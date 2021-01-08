%% This script generate the relative error curve for each of the method to 
%% show the method is convergent as T goes by

clc;
clear;
clf;

% Parameter values:
K = 50;
S_max = 3*K;
r = 0.01;
sigma = 0.25;
dt = 3/48000;
ds = S_max/500;
T = 0:dt:3.0;

% % check for stability of explicit method
% dt = 3/48000;
% ds = S_max/500;
% S_max2 = S_max^2;
% sigma2 = sigma^2;
% dt_ds2 = dt./ds.^2;
% Constriant = 1/(sigma2*S_max2);
% 
% % Verify the stability condition satisfied for explicit method
% dt_ds2 <= Constriant

%% Explicit

% Solution for the Explicit
[S,V_explicit,tao] = Explicit_B_S(K,S_max,r,T(end),sigma,ds,dt,'CALL');
% Solution for the Analytic formula
V_exact = Exact_B_S(K,S_max,r,T,sigma,ds,'CALL');

% Reverse the order of V_exact
V_exact = fliplr(V_exact);

% Initilzed the error vector
L_explicit_error = zeros(length(T),1);

% Calculate the relative error over time t from 0 to 1
for i  = 1:length(T)
    L_explicit_error(i) = vecnorm(V_exact(:,i)-V_explicit(:,i),1) / vecnorm(V_exact(:,i),1);
end

%% Semi-implicit

% Solution for the Semi-Implicit
[S,V_semi,A_semi] = Semi_implicit_B_S(K,S_max,r,T(end),sigma,ds,dt,'CALL');

% check stability of method (stable if infiniry norm <= 1)
Ainv_semi = inv(A_semi);
normi_semi = norm(Ainv_semi,inf);

% Initilzed the error vector
L_semi_error = zeros(length(T),1);

% Calculate the relative error over time t from 0 to 1
for i  = 1:length(T)
    L_semi_error(i) = vecnorm(V_exact(:,i)-V_semi(:,i),1) / vecnorm(V_exact(:,i),1);
end

%plot the time versus relative error curve
figure(1)
plot(T,L_semi_error,'-b')
xlabel('time (t)')
ylabel('Relative Error in L-1 Norm')
legend('Relative error for \Delta t = 6.25e-05 and \Delta x = 0.3')


%% Fully-implicit

% Solution for the Fully-Implicit
[S,V_fully,A_fully] = Fully_Implicit_B_S(K,S_max,r,T(end),sigma,ds,dt,'CALL');

% check stability of method (stable if infiniry norm <= 1)
Ainv_fully = inv(A_fully);
normi_fully = norm(Ainv_fully,inf);

% Initilzed the error vector
L_fully_error = zeros(length(T),1);

% Calculate the relative error over time t from 0 to 1
for i  = 1:length(T)
    L_fully_error(i) = vecnorm(V_exact(:,i)-V_fully(:,i),1) / vecnorm(V_exact(:,i),1);
end

%plot the time versus relative error curve
figure(2)
plot(T,L_fully_error,'-b')
xlabel('time (t)')
ylabel('Relative Error in L-1 Norm')
legend('Relative error for \Delta t = 6.25e-05 and \Delta x = 0.3')


%% CN

% Solution for the CN
[S,V_CN,A_CN,B_CN] = CN_B_S(K,S_max,r,T(end),sigma,ds,dt,'CALL');

% check stability of method (stable if infiniry norm <= 1)
Ainv_CN = inv(A_CN);
normi_CN = norm(Ainv_CN.*B_CN,inf);

% Initilzed the error vector
L_CN_error = zeros(length(T),1);

% Calculate the relative error over time t from 0 to 1
for i  = 1:length(T)
    L_CN_error(i) = vecnorm(V_exact(:,i)-V_CN(:,i),1) / vecnorm(V_exact(:,i),1);
end

%plot the time versus relative error curve
figure(3)
plot(T,L_CN_error,'-b')
xlabel('time (t)')
ylabel('Relative Error in L-1 Norm')
legend('Relative error for \Delta t = 6.25e-05 and \Delta x = 0.3')
