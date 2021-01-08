%% Accuracy Table and Accuracy plot

clc;
clear;
clf;

% Example: 

% S = 10 20 30 ... 90;
Stock = 10:5:90;
% Strike price K 
K = 50;
% S_max = 3*k; as Willmot advised
S_max = 3*K;
% risk free interest rate
r = 0.05;
% T = 1 years
T = 3;
% sigma = 0.20
sigma = 0.25;
% Asset step M = 500
M = 500;
% take time step = 50000 so that the explicit method will not be un stable
N = 50000;

% check for stability of explicit method
dt = T/N;
ds = S_max/M;
S_max2 = S_max^2;
sigma2 = sigma^2;
dt_ds2 = dt./ds.^2;
Constriant = 1/(sigma2*S_max2);

% Verify the stability condition satisfied for explicit method
dt_ds2 <= Constriant

%% Start

%% #1 analytic solution for each stock price
% Solution for the Analytic formula 0.031803 seconds)
tic,result_exact = zeros(length(Stock),1);

for i = 1:length(Stock)
    V_exact = Exact_B_S(K,S_max,r,T,sigma,ds,'PUT');
    % Reverse the order of V_exact
    V_exact = fliplr(V_exact);
    result_exact(i) = interp1(0:ds:S_max, V_exact(:,1),Stock(i));
end,toc

%% #2 Explicit solution for each stock price

% Solution for the Explicit (0.396048 seconds.)
tic,result_explicit = zeros(length(Stock),1);
[S,V_explicit,tao] = Explicit_B_S(K,S_max,r,T,sigma,ds,dt,'PUT');
for i = 1:length(Stock)
    result_explicit(i) = interp1(S, V_explicit(:,1),Stock(i));
end,toc

%% #3 Semi-implicit solution for each stock price

% Solution for the Semi-implicit (18.151862 seconds.)
tic,result_semi_explicit = zeros(length(Stock),1);
[S,V_semi_explicit,tao] = Semi_implicit_B_S(K,S_max,r,T,sigma,ds,dt,'PUT');
for i = 1:length(Stock)
    result_semi_explicit(i) = interp1(S, V_semi_explicit(:,1),Stock(i));
end,toc

%% #4 Fully-implicit solution for each stock price
% Solution for the Fully-implicit (16.484444 seconds.)
tic,result_fully_explicit = zeros(length(Stock),1);
[S,V_fully_implicit,tao] = Fully_Implicit_B_S(K,S_max,r,T,sigma,ds,dt,'PUT');
for i = 1:length(Stock)
    result_fully_explicit(i) = interp1(S, V_fully_implicit(:,1),Stock(i));
end,toc

%% #5 Crank-nicolson solution for each stock price
% Solution for the Crank-nicolson (19.718434 seconds.)
tic,result_CN = zeros(length(Stock),1);
[S,V_CN,tao] = CN_B_S(K,S_max,r,T,sigma,ds,dt,'PUT');
for i = 1:length(Stock)
    result_CN(i) = interp1(S, V_CN(:,1),Stock(i));
end,toc

%% #6 Monte_carlo solution for each stock price
% Solution for the Monte-Carlo (3.417380 seconds.)
path = 1e7;

tic,result_MC = zeros(length(Stock),1);
for i = 1:length(Stock)
    result_MC(i) = Monte_carlo(Stock(i),K,r,sigma,T,path);
end,toc

%% Finally, Calculated the relative error for each stock price

%% #1 Relative Error for Explicit solution for each stock price

% Initilzed the error vector
Absolute_error_explicit = zeros(length(Stock),1);

% Calculate the relative error over time t from 0 to 1
for i  = 1:length(Stock)
    Absolute_error_explicit(i) = abs(result_exact(i) - result_explicit(i));
end

%% #2 Relative Error for Semi-implicit solution for each stock price

% Initilzed the error vector
Absolute_error_semi_implicit = zeros(length(Stock),1);

% Calculate the relative error over time t from 0 to 1
for i  = 1:length(Stock)
    Absolute_error_semi_implicit(i) = abs(result_exact(i) - result_semi_explicit(i));
end

%% #3 Relative Error for Fully-implicit solution for each stock price
% Initilzed the error vector
Absolute_error_Fully_implicit = zeros(length(Stock),1);

% Calculate the relative error over time t from 0 to 1
for i  = 1:length(Stock)
    Absolute_error_Fully_implicit(i) = abs(result_exact(i) - result_fully_explicit(i));
end

%% #4 Relative Error for Crank-Nicolson solution for each stock price

% Initilzed the error vector
Absolute_error_CN = zeros(length(Stock),1);

% Calculate the relative error over time t from 0 to 1
for i  = 1:length(Stock)
    Absolute_error_CN(i) = abs(result_exact(i) - result_CN(i));
end

%% #4 Relative Error for Monte_carlo solution for each stock price

% Initilzed the error vector
Absolute_error_MC = zeros(length(Stock),1);

% Calculate the relative error over time t from 0 to 1
for i  = 1:length(Stock)
    Absolute_error_MC(i) = abs(result_exact(i) - result_MC(i));
end


% All together:

figure(4)
hold on
plot(Stock,Absolute_error_explicit,'-r')
plot(Stock,Absolute_error_semi_implicit,'-b')
plot(Stock,Absolute_error_Fully_implicit,'-p')
plot(Stock,Absolute_error_CN,'-y')
plot(Stock,Absolute_error_MC,'g')
xlabel('Stock Price')
ylabel('Absolute Error')
legend('explicit', 'semi_implicit','Fully_implicit','Crank Nicolson', 'Monte Carlo')
hold off