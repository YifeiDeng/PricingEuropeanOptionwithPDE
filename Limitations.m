%% Some Results for Executes summary presentation

%% Limitation for explicit method

%% #1

% the right one
[S_1,V_correct,tao] = Explicit_B_S(60,100,0.05,1,0.2,1,0.001,'CALL');
subplot(2,1,1)
plot(S_1,V_correct(:,1),'r')
xlabel 'Stock Price'
ylabel 'Option Price'
title('The Stable Solution for Explicit Algorithm')

% the incorrect one
subplot(2,1,2)
[S_2,V_incorrect,tao] = Explicit_B_S(60,100,0.05,1,0.2,1,0.01,'CALL');
plot(S_2,V_incorrect(:,1),'r')
xlabel 'Stock Price'
ylabel 'Option Price'
title('The Unstable Solution for Explicit Algorithm')

%% #2 limitation of the other algorithms

% limitation of the other algorithm

% Parameter value from the paper:
K = 100;
S_max = 130;
r = 0.12;
sigma = 0.10;
dt = 0.02;
ds = 2;
T = 0:dt:1.0;

% % check for stability of explicit method
% dt = T/N;
% ds = S_max/M;
% S_max2 = S_max^2;
% sigma2 = sigma^2;
% dt_ds2 = dt./ds.^2;
% Constriant = 1/(sigma2*S_max2);
% 
% % Verify the stability condition satisfied for explicit method
% dt_ds2 <= Constriant

tic,[S,V,tao] = CN_B_S(K,S_max,r,T(end),sigma,ds,dt,'CALL');toc;
% Solution for the Analytic formula
V_exact = Exact_B_S(K,S_max,r,T,sigma,ds,'CALL');

% Initilzed the error vector
L_1_error_1 = zeros(length(T),1);

% Reverse the order of V_exact
V_exact = fliplr(V_exact);

% Calculate the relative error over time t from 0 to 1
for i  = 1:length(T)
    L_1_error_1(i) = vecnorm(V_exact(:,i)-V(:,i),1) / vecnorm(V_exact(:,i),1)*100;
end

L_1_error_1(end) = 0;

figure(2)
plot(T,L_1_error_1,'-b')
xlabel('time (t)')
ylabel('Percentage Error')
legend('Percentage Error for N = 1600000 and M = 1200')