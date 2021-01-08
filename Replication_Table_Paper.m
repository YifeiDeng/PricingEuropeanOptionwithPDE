%% Replication of tables in the paper

%% ===================clear environment====================================

clear
close all
clc   
clf

%% Table 1

% stock price vector
Stock = [4 8 10 16 20];

% Other parameter consider
S_min = 0;
S_max = 20;
T = 0.25;
K = 10;
r = 0.1;
sigma = 0.4;
N = 41000;
M = 1000;
ds = S_max / M;
dt = T / N;

% Explicit method 
tic,[S,V,tao] = Explicit_B_S(K,S_max,r,T,sigma,ds,dt,'CALL');,toc
stock_Ex_1 = interp1(S, V(:,1),Stock(1));
stock_Ex_2 = interp1(S, V(:,1),Stock(2));
stock_Ex_3 = interp1(S, V(:,1),Stock(3));
stock_Ex_4 = interp1(S, V(:,1),Stock(4));
stock_Ex_5 = interp1(S, V(:,1),Stock(5));

% Exact formula
tic,V_exact = Exact_B_S(K,S_max,r,T,sigma,ds,'CALL');,toc
stock_E_1 = interp1(S, V_exact(:,1),Stock(1));
stock_E_2 = interp1(S, V_exact(:,1),Stock(2));
stock_E_3 = interp1(S, V_exact(:,1),Stock(3));
stock_E_4 = interp1(S, V_exact(:,1),Stock(4));
stock_E_5 = interp1(S, V_exact(:,1),Stock(5));

% Semi-implict method
tic,[S,V_semi,tao] = Semi_implicit_B_S(K,S_max,r,T,sigma,ds,dt, 'CALL');,toc
stock_s_1 = interp1(S, V_semi(:,1),Stock(1));
stock_s_2 = interp1(S, V_semi(:,1),Stock(2));
stock_s_3 = interp1(S, V_semi(:,1),Stock(3));
stock_s_4 = interp1(S, V_semi(:,1),Stock(4));
stock_s_5 = interp1(S, V_semi(:,1),Stock(5));

% fully implict method
tic,[S,V_im,tao] = Fully_Implicit_B_S(K,S_max,r,T,sigma,ds,dt, 'CALL');,toc
stock_i_1 = interp1(S, V_im(:,1),Stock(1));
stock_i_2 = interp1(S, V_im(:,1),Stock(2));
stock_i_3 = interp1(S, V_im(:,1),Stock(3));
stock_i_4 = interp1(S, V_im(:,1),Stock(4));
stock_i_5 = interp1(S, V_im(:,1),Stock(5));

% CN method
tic,[S,V_CN,A,B] = CN_B_S(K,S_max,r,T,sigma,ds,dt,'CALL');,toc
stock_cn_1 = interp1(S, V_CN(:,1),Stock(1));
stock_cn_2 = interp1(S, V_CN(:,1),Stock(2));
stock_cn_3 = interp1(S, V_CN(:,1),Stock(3));
stock_cn_4 = interp1(S, V_CN(:,1),Stock(4));
stock_cn_5 = interp1(S, V_CN(:,1),Stock(5));

%% Table 2

% stock price vector
Stock = [4 8 10 16 20];

% Other parameter consider
S_min = 0;
S_max = 20;
T = 0.25;
K = 10;
r = 0.1;
sigma = 0.4;
N = 41000;
M = 1000;
ds = S_max / M;
dt = T / N;

% Explicit method 
tic,[S,V,tao] = Explicit_B_S(K,S_max,r,T,sigma,ds,dt,'CALL');,toc
stock_Ex_1 = interp1(S, V(:,1),Stock(1));
stock_Ex_2 = interp1(S, V(:,1),Stock(2));
stock_Ex_3 = interp1(S, V(:,1),Stock(3));
stock_Ex_4 = interp1(S, V(:,1),Stock(4));
stock_Ex_5 = interp1(S, V(:,1),Stock(5));

% Exact formula
tic,V_exact = Exact_B_S(K,S_max,r,T,sigma,ds,'CALL');,toc
stock_E_1 = interp1(S, V_exact(:,1),Stock(1));
stock_E_2 = interp1(S, V_exact(:,1),Stock(2));
stock_E_3 = interp1(S, V_exact(:,1),Stock(3));
stock_E_4 = interp1(S, V_exact(:,1),Stock(4));
stock_E_5 = interp1(S, V_exact(:,1),Stock(5));