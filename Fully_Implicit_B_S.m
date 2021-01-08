% Fully Implicit difference scheme for a European-type option in Black-Scholes
% equations
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
% A: matrix A in the matrix notation of the algorithm for convergence
% analysis
%--------------------------------------------------------------------------
function [S,V,A] = Fully_Implicit_B_S(K,S_max,r,T,sigma,ds,dt,option_type)
%% Grid setup and Boundary condition %%
 
% Price steps
M = round((S_max) / ds);
ds = (S_max) / M;
% Time steps
N = round(T / dt);
dt = T / N;
% Initialized the matrix for option value
V = zeros(M+1,N+1);
% Discretize the asset S
S = linspace(0,S_max,M+1)';
% Discretize the time t
t = 0:dt:T;
 
% Specify the type of the vallina option
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
        
        % Impose boundary condition under the European call: 
        % V(0,t) = K*exp(-rt) & V(S,t) = 0 as S goes to infinity
        V(1, :) = K * exp(-r*t(end:-1:1));
        V(end, :) = 0;
end
 
% Set coefficient
veci = 0:M;
 
a = 0.5 * (r * dt * veci - sigma^2 * dt * (veci .^2));
b = 1 + sigma^2 * dt * (veci .^2) + r * dt;
c = -0.5 * (r * dt * veci + sigma^2 * dt * (veci.^2)); 
 
% Construction of solution matrix
A = diag(a(3:M),-1) + diag(b(2:M)) + diag(c(2:M-1),1);
 
% Apply LU-decomposition
[L,U] = lu(A);
 
% Construct the vector for remaining coefficient
z = zeros(size(A,2),1);
 
% Perform the semi-implicit scheme (time: n = 1...N, i = 2...m for V(i,n))
for n = N:-1:1
    z(1) = - a(2) * V(1,n);
    z(end) = - c(end) * V(end,n);   
    V(2:M,n) =  U \ (L \ (V(2:M,n+1) + z));
end

end