% Function to perform LASSO regression using Alternating Direction Method
% of Multipliers.
%
% arg min_{B} 0.5*||X - A*B||_{2}^{2} + alpha*||B||_{1}
%
% Usage:- [B,cost] = lasso_admm(X, A, alpha)
%
% where:-
%         b = bias vector
%         lambda = weighting on the l1 penalty
%
%         x = solution
%
% Written by Simon Lucey 2012

function [x,cost,c_n] = lasso_admm_lap(y, A, L, AA, LL, alpha)

%keyboard;

% Get dimensions of B
c = size(y,2);
r = size(A,2); 

lagragian = randn(r,c); % Initialize Lagragian to be nothing (seems to work well)
rho = 1e-3; % Set rho to be quite low to start with
maxIter = 20; % Set the maximum number of iterations (make really big to ensure convergence)
I = speye(r); % Set the sparse identity matrix
maxRho = 5; % Set the maximum mu
b = randn(r,c); % Initialize C randomly

% Set the fast soft thresholding function
fast_sthresh = @(x,th) sign(x).*max(abs(x) - th,0);

% Set the norm functions
norm2 = @(x) x(:)'*x(:);
norm1 = @(x) sum(abs(x(:))); 

% AA = A'*A;
% LL = L'*L;
A1 = A'*y;

% cost = [];
cost = 0;
for n = 1:maxIter

    % Solve sub-iproblem to solve B
    x = (AA+rho*LL)\(A1 + L'*(rho*b - lagragian) ); 

    % Solve sub-problem to solve C
    b = fast_sthresh(L*x + lagragian/rho, alpha/rho); 

    % Update the Lagrangian
    lagragian = lagragian + rho*(L*x - b);  

    %pause; 

    % Section 3.3 in Boyd's book describes strategies for adapting rho
    % main strategy should be to ensure that
    rho = min(maxRho, rho*1.1); 

    % get the current cost
    cost(n) = 0.5*norm2(y - A*x) + alpha*norm1(L*x);

%     figure(1)
%     subplot(311); plot(x);
%     subplot(312); plot(b);
%     subplot(313); plot(lagragian);
%     pause(0.01);
    
    %if n==10 keyboard; end
end
c_n=0;
% if ~isempty(find(abs(diff(cost))<0.01*abs(cost(2)-cost(1)),1))
%     c_n = find(abs(diff(cost))<0.01*abs(cost(2)-cost(1)),1)+1;
% else
%     c_n = [];
% end

