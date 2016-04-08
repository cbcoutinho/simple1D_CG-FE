function value = basis_1D(order, N, x, dx)
%BASIS_1D Calculates the value of a 1D isoparametric basis function
%   Calculates the value of a 1D isoparametric basis function. An
%   isoparametric basis function has a fixed range of [-1, 1], which is
%   required for numerical integration using Guassian-Legendre quadrature.
% 
%   The variables are:
%       order:  Order of the basis function (1 = linear, 2 = quadratic, etc.)
%       N:      Basis function number
%       x:      Coordinate at which to calculate the basis function
%       dx:     Degree of differentiation (0 = none, 1 = d/dx, 2 = d2/dx2, etc.)
%       value:  The output value of the basis function

switch order
    case 1  % Linear
        if N > 2
            error('Point requested does not fall in available basis functions')
        else
            value = linear(N, x, dx);
        end
    case 2 % Quadratic
        if N > 3
            error('Point requested does not fall in available basis functions')
        else
            value = quadratic(N, x, dx);
        end
    case 3 % Cubic
        if N > 4
            error('Point requested does not fall in available basis functions')
        else
            value = cubic(N, x, dx);
        end
    otherwise
        error(['order = ', num2str(order), ' selected that is not supported'])
end

end

function value = linear(N, x, dx)
%LINEAR Calculates the value of the linear isoparametric basis function
% 
%   The variables are:
%       N:      The basis number
%       x:      The value at which to calculate the basis function
%       value:  The output value of the basis function

switch N
    case 1
        if dx == 0
            value = 0.5 - 0.5*x;
        else
            value = -0.5*ones(length(x),1);
        end
    case 2
        if dx == 0
            value = 0.5 + 0.5*x;
        else
            value = 0.5*ones(length(x),1);
        end
end

end

function value = quadratic(N, x, dx)
%QUADRATIC Calculates the value of the quadratic isoparametric basis
%   function
% 
%   The variables are:
%       N:      The basis number
%       x:      The value at which to calculate the basis function
%       value:  The output value of the basis function

switch N
    case 1
        if dx == 0
            value = -0.5*x + 0.5*x.^2;
        else
            value = -0.5 + x;
        end
    case 2
        if dx == 0
            value = 1 - x.^2;
        else
            value = -2*x;
        end
    case 3
        if dx == 0
            value = 0.5*x + 0.5*x.^2;
        else
            value = 0.5 + x;
        end
end

end

function value = cubic(N, x, dx)
%CUBIC Calculates the value of the cubic isoparametric basis function
% 
%   The variables are:
%       N:      The basis number
%       x:      The value at which to calculate the basis function
%       value:  The output value of the basis function

switch N
    case 1
        if dx == 0
            value = -0.0625 + 0.0625*x + 0.5625*x.^2 - 0.5625*x.^3;
        else
            value = 0.0625 + 2*0.5625*x - 3*0.5625*x.^2;
        end
    case 2
        if dx == 0
            value = 0.5625 - 1.6875*x - 0.5625*x.^2 + 1.6875*x.^3;
        else
            value = -1.6875 - 2*0.5625*x + 3*1.6875*x.^2;
        end
    case 3
        if dx == 0
            value = 0.5625 + 1.6875*x - 0.5625*x.^2 - 1.6875*x.^3;
        else
            value = 1.6875 - 2*0.5625*x - 3*1.6875*x.^2;
        end
    case 4
        if dx == 0
            value = -0.0625 - 0.0625*x + 0.5625*x.^2 + 0.5625*x.^3;
        else
            value = -0.0625 + 2*0.5625*x + 3*0.5625*x.^2;
        end
end

end