% This is wrapper script is used to build the required matricies and
% vectors for a quick FEM analysis. The function attempting to be solved is
% as follows:
%       a*du/dx = cos(x) on [0, 2*pi]
%
%   where
%       a:  Scalar parameter
%       u:  Scalar value of the function
%       x:  Independent coordinate
%
%   The theoretical solution to this ODE is:
%       u = sin(x)/a

L   = [0, 2*pi];    % Domain bounds
u0  = [0, 0];       % Boundary Conditions
numElem = 7;        % Number of Elements (Currently default linear)
a = 1;              % Parameter value

coords = linspace(L(1), L(2), numElem+1)';
coordPairs = zeros(numElem,2);

for ii = 1:numElem
    coordPairs(ii,:) = coords(ii:ii+1);
end

A = zeros(numElem+1);
b = zeros(numElem+1,1);

for kk = 1:numElem
    Ie  = zeros(2);
    X   = @(x)basis_1D(1,1,x,0)*coordPairs(kk,1) + ...
                basis_1D(1,2,x,0)*coordPairs(kk,2); % X = x(s)
    J   = @(x)basis_1D(1,1,x,1)*coordPairs(kk,1) + ...
                basis_1D(1,2,x,1)*coordPairs(kk,2); % J = dx/ds
    for ii = 1:2
        for jj = 1:2
            fun1    = @(s) basis_1D(1,ii,s,1)./J(s);
            fun2    = @(s) basis_1D(1,jj,s,0);
            fun     = @(s) fun1(s)*fun2(s).*J(s);
            
            Ie(ii,jj) = -a*lgquad(fun,[],-1,1);
            
            II = (kk-1)+ii;
            JJ = (kk-1)+jj;
            A(II,JJ) = A(II,JJ) + Ie(ii,jj);
        end
        fun1    = @(s) cos(X(s));
        fun2    = @(s) basis_1D(1,ii,s,0);
        fun     = @(s) fun1(s)*fun2(s).*J(s);
        b((kk-1)+ii) = b((kk-1)+ii) + lgquad(fun,[],-1,1);
    end
end

A = A(2:end-1,2:end-1);
b = b(2:end-1);
x = A\b;

u = [u0(1);x;u0(2)];

figure(1); clf; hold on; grid on;

% Plot exact solution
plot(linspace(L(1),L(2),100),sin(linspace(L(1),L(2),100))/a,'-')

% Plot the solution of each element separately
for ii = 1:numElem
    plot(coordPairs(ii,:),u(ii:ii+1),'rs-')
end

ax = gca;
ax.XTick = [0 pi/2 pi 3*pi/2 2*pi];
ax.XTickLabel = {'0','\pi/2','pi','3\pi/2','2\pi'};
legend('Exact Solution','FEM Solution')


