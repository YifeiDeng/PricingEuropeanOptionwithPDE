% This is the program using Monte-Carlo simlulation for pricng a
% European type of put opiton
%--------------------------------------------------------------------------
% INPUTS:
%
%   S:      Stock price
%   K:      strike price
%   r:      risk-free interest rate
%   sigma:  volatility
%   T:   time to maturity
%   M: number of path to be specified in simulation
%--------------------------------------------------------------------------
%
% OUTPUT:
%
% Put_price: price of a European put option
%
%--------------------------------------------------------------------------
function Put_price = Monte_carlo(S,K,r,sigma,T,M)
% the stock price at expiry T
S_T = S*exp((r-0.5*sigma^2)*T + sigma*sqrt(T)*randn(M,1));

% Evaluate the Put option options
option_values = max(K-S_T,0);

% Discount Value under assumption
present_vals = exp(-r*T)*option_values;
% Take the average
put_value = mean(present_vals);

% Return the European put price of the option
Put_price = put_value;
end
