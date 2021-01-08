% Explicit difference scheme for a European type option in Black-Scholes
% equation
%--------------------------------------------------------------------------
% INPUTS:
%
%   K:      strike price
%   S_max:  the largest stock price
%   r:      risk-free interest rate
%   T:      time to maturity
%   sigma:  volatility
%   ds:     size for price step
%   dt:     size for time step
%   option_type:  option types, type either'CALL' or 'PUT'
%--------------------------------------------------------------------------
%
% OUTPUT:
%
% S: discretized price domain
% V: option value at each time and each step
% tao: discretized time series
%--------------------------------------------------------------------------
function [S,V,tao] = Explicit_B_S(K,S_max,r,T,sigma,ds,dt,option_type)
%% Grid setup and Boundary condition %%
 
% Price steps
M = round((S_max) / ds);
ds = (S_max) / M;
% Time steps
N = round(T / dt);
dt = T / N;
% Initialized the matrix for option value
V(1:M+1,1:N+1) = 0.0;
% Discretize the asset S
S = linspace(0,S_max,M+1);
tao = (0:N) * dt;
% Discretize the asset t
t = 0:dt:T;
 
switch option_type
    case 'CALL'
        % Initial condition of European call Payoff at expiry T: V(S,T) =
        % (S-K)+
        V(:, end) = max(S - K,0);
 
        % Impose boundary condition under the European call: 
        % V(0,t) = 0 & V(S,t) = S - K*exp(-r(T-t)) as S goes to infinity
        V(1, :) = 0;
        V(end, :) = S_max - K * exp(-r*t(end:-1:1));
    case 'PUT'
        % Initial condition of European put Payoff at expiry T: V(S,T) =
        % (K-S)+
        V(:, end) = max(K - S,0);
        
        % Impose boundary condition under the European put: 
        % V(0,t) = K*exp(-rt) & V(S,t) = 0 as S goes to infinity
        V(1, :) = K * exp(-r*t(end:-1:1));
        V(end, :) = 0;
end
 
% Set coefficient
veci = 0:M;
a = 0.5 * dt * (sigma^2 * veci - r).*veci;
b = 1 - dt * (sigma^2 * veci.^2 + r);
c = 0.5 * dt * (sigma^2 * veci + r).*veci;
 
% Perform the explicit scheme (time: n = 1...N, i = 2...m for V(i,n))
for n = N:-1:1
    for i = 2:M
        V(i,n) = a(i) * V(i-1,n+1)...
                + b(i) * V(i,n+1)...
                + c(i) * V(i+1,n+1);
    end
end
 
end