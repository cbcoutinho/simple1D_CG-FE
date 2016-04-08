function out = lgquad( fun, N, a, b )
%LGQUAD Computes the integral of a function
%   Guass-Legendre quadrature is used to compute the integral of a
%   continuous polynomial function. The function first calls a sub-function
%   to calculate the coordinates and weights required in the integral, and
%   then calculates the integral numerically.
%
%   If N is not provided (empty), then a wrapper is used to incrementally
%   increase the order of integration until convergence
%
%   The variables are:
%       fun:    A function object which accepts only a single input f(x)
%       N:      Number of Legendre points used in the integral
%       a:      Lower limit of integration
%       b:      Upper limit of integration
%       out:    Calculated integral of fun on domain [a, b]

if ~isempty(N)
    % If N is provided, then calculate normally.
    out = lgquadCalc(fun, N, a, b);
else
    % If N is not provided, wrap around lgquadCalc and incrementally
    % increase the order of integration until convergence.
    epsilon = sqrt(eps);
    diff = 1;
    
    N = 1;
    while diff > epsilon
        
        out = lgquadCalc(fun, N, a, b);
        
        if N == 1
            out_old = out;
        else
            diff = abs(out-out_old);
            out_old = out;
        end
        
        N = N+1;
    end
end

end

function out = lgquadCalc(fun, N, a, b)

[x,c] = lgwt(N,-1,1);
len = length(x);

out = 0;

x_star = (b-a)/2.*x+(b+a)/2;
c_star = (b-a)/2.*c;

for ii = 1:len
    out = out + c_star(ii) * fun(x_star(ii));
end

end