% Exact solution for a European type option in Black-Scholes
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
%   option_type:  option types, type either'CALL' or 'PUT'
%--------------------------------------------------------------------------
%
% OUTPUT:
%
% V: option value at each time and each step
%--------------------------------------------------------------------------
function V = Exact_B_S(K,S_max,r,T,sigma,ds,option_type)
%% Grid setup and Boundary condition %%
% T
times = length(T);

% Price steps
M = round((S_max) / ds);
ds = (S_max) / M;

% Price grids
S = linspace(0,S_max,M+1);

% Initialized the option value
V = zeros(length(S),length(T));

switch option_type
    case 'CALL'
        % Exact European call solution using classical black scholes formula:
        for i = 1:times
            d1 = (log(S ./ K) + (r + 1/2 * sigma^2) * (T(i))) / (sigma * sqrt(T(i)));
            d2 = (log(S ./ K) + (r - 1/2 * sigma^2) * (T(i))) / (sigma * sqrt(T(i)));
            V(:,i) = S .* normcdf(d1) - K .* exp(-r * (T(i))) * normcdf(d2);
        end
    case 'PUT'
        % Exact European put solution using classical black scholes formula:
        for i = 1:times
            d1 = (log(S ./ K) + (r + 1/2 * sigma^2) * (T(i))) / (sigma * sqrt(T(i)));
            d2 = (log(S ./ K) + (r - 1/2 * sigma^2) * (T(i))) / (sigma * sqrt(T(i)));
            V(:,i) = K .* exp(-r * (T(i))) * normcdf(-d2) - normcdf(-d1) .* S;
        end  
end


end

