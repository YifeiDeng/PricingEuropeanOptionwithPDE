% Crank-Nicolson scheme for a European type option in Black-Scholes
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
% B: matrix B in the matrix notation of the algorithm for convergence
% analysis
%--------------------------------------------------------------------------
function [S,V, A,B] = CN_B_S(K,S_max,r,T,sigma,ds,dt,option_type)
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
S = 0:ds:S_max;
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
        
        % Impose boundary condition under the European Put: 
        % V(0,t) = K*exp(-rt) & V(S,t) = 0 as S goes to infinity
        V(1, :) = K * exp(-r*t(end:-1:1));
        V(end, :) = 0;
end

% Set coefficients
veci = 0:M;
sigma2 = sigma*sigma;

ai = (dt/4)*(sigma2*(veci.^2) - r*veci);
bi = -(dt/2)*(sigma2*(veci.^2) + r);
ci = (dt/4)*(sigma2*(veci.^2) + r*veci);

% Form the tridiagonal matrix
A = -diag(ai(3:M),-1) + diag(1-bi(2:M)) - diag(ci(2:M-1),1);
B = diag(ai(3:M),-1) + diag(1+bi(2:M)) + diag(ci(2:M-1),1);

% Apply LU-decomposition
[L,U] = lu(A);

% Solve at each node
aux = zeros(size(B,2),1);

for n = N:-1:1
    if length(aux)==1
        aux = ai(2)*(V(1,n)+V(1,n+1)) + ...
              ci(end)*(V(end,n)+V(end,n+1));
    else
        aux(1) = ai(2)*(V(1,n)+V(1,n+1));
        aux(end) = ci(end)*(V(end,n)+V(end,n+1));
    end
    V(2:M,n) = U\(L\(B*V(2:M,n+1) + aux));
end


end

