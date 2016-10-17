function [x_hat, lambda_opt] = dsvd (A, y, lambda, SNR, compute_params)
% DSVD reconstruction method.
%
% This routine, by victor.suarez.gutierrez@urjc.es
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUT:
% A: transfer matrix.
% y: torso potentials.
% lambda: regularization parameter.
% SNR: signal noise rate (dB).
% compute_params: True whether reg_params need to be calculated. Otherwise False. 

%OUTPUT:
% x_hat: estimated epicardial potentials.
% lambda_opt: optimal regularization parameters (one if global or an array by instant).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Calculation:
[U,sm,V]   = csvd(A);
   
if ~compute_params
    lambda_opt = lambda;
    reg_corner = sqrt(lambda);
    for t=1:length(y(1,:))
        fprintf('DSVD Instant: t = %i.\n', t);
        x_hat(:,t) = dgsvd (U, sm, V, y(:,t), reg_corner(t));
    end
else
    for t=1:length(y(1,:))
        fprintf('DSVD Instant: t = %i.\n', t);
        [reg_corner, rho, eta, reg_param] = l_curve (U, sm, sm, y(:,t), 'dsvd', SNR);
        lambda_opt(1,t) = reg_corner^2;
        x_hat(:,t) = dgsvd (U, sm, V, y(:,t), reg_corner);
    end
end   
end


function [U,sm,X,V,W] = cgsvd(A,L)
%CGSVD Compact generalized SVD (GSVD) of a matrix pair in regularization problems.
%
% sm = cgsvd(A,L)
% [U,sm,X,V] = cgsvd(A,L) ,  sm = [sigma,mu]
% [U,sm,X,V,W] = cgsvd(A,L) ,  sm = [sigma,mu]
%
% Computes the generalized SVD of the matrix pair (A,L). The dimensions of
% A and L must be such that [A;L] does not have fewer rows than columns.
%
% If m >= n >= p then the GSVD has the form:
%    [ A ] = [ U  0 ]*[ diag(sigma)      0    ]*inv(X)
%    [ L ]   [ 0  V ] [      0       eye(n-p) ]
%                     [  diag(mu)        0    ]
% where
%    U  is  m-by-n ,    sigma  is  p-by-1
%    V  is  p-by-p ,    mu     is  p-by-1
%    X  is  n-by-n .
%
% Otherwise the GSVD has a more complicated form (see manual for details).
%
% A possible fifth output argument returns W = inv(X).
 
% Reference: C. F. Van Loan, "Computing the CS and the generalized 
% singular value decomposition", Numer. Math. 46 (1985), 479-491. 
 
% Per Christian Hansen, DTU Compute, August 22, 2009. 
 
% Initialization.
[m,n] = size(A); [p,n1] = size(L);
if (n1 ~= n)
  error('No. columns in A and L must be the same')
end
if (m+p < n)
  error('Dimensions must satisfy m+p >= n')
end

% Call Matlab's GSVD routine.
[U,V,W,C,S] = gsvd(full(A),full(L),0);

if (m >= n)
  % The overdetermined or square case.
  %sm = [diag(C(1:p,1:p)),diag(S(1:p,1:p))];
  q = min(p,n);
  sm = [diag(C(1:q,1:q)),diag(S(1:q,1:q))];
  if (nargout < 2) 
    U = sm; 
  else 
    % Full decomposition. 
    X = inv(W'); 
  end
else
  % The underdetermined case.
  sm = [diag(C(1:m+p-n,n-m+1:p)),diag(S(n-m+1:p,n-m+1:p))]; 
  if (nargout < 2) 
    U = sm; 
  else 
    % Full decomposition. 
    X = inv(W');
    X = X(:,n-m+1:n); 
  end
end

if (nargout==5), W = W'; end
end


function [x_lambda,rho,eta] = dgsvd (U,s,V,b,lambda)
%DSVD Damped SVD and GSVD regularization.
%
% [x_lambda,rho,eta] = dsvd(U,s,V,b,lambda)
% [x_lambda,rho,eta] = dsvd(U,sm,X,b,lambda) ,  sm = [sigma,mu]
%
% Computes the damped SVD solution defined as
%    x_lambda = V*inv(diag(s + lambda))*U'*b .
% If lambda is a vector, then x_lambda is a matrix such that
%    x_lambda = [ x_lambda(1), x_lambda(2), ... ] .
% U, s, and V must be computed by the csvd function.
%
% If sm and X are specified, then the damped GSVD solution is computed:
%    x_lambda = X*[ inv(diag(sigma + lambda*mu)) 0 ]*U'*b
%                 [            0                 I ]
% U, sm, and X must be computed by the cgsvd function.
%
% The solution norm (standard-form case) or seminorm (general-form
% case) and the residual norm are returned in eta and rho.

% Reference: M. P. Ekstrom & R. L. Rhoads, "On the application of
% eigenvector expansions to numerical deconvolution", J. Comp.
% Phys. 14 (1974), 319-340.
% The extension to GSVD is by P. C. Hansen.

% Per Christian Hansen, DTU Compute, April 14, 2003.

% Initialization.
if (min(lambda)<0)
  error('Illegal regularization parameter lambda')
end
m = size(U,1);
n = size(V,1);
[p,ps] = size(s);
beta = U(:,1:p)'*b;
ll = length(lambda); x_lambda = zeros(n,ll);
rho = zeros(ll,1); eta = zeros(ll,1);

% Treat each lambda separately.
if (ps==1)

  % The standard-form case.
  for i=1:ll
    x_lambda(:,i) = V(:,1:p)*(beta./(s + lambda(i)));
    rho(i) = lambda(i)*norm(beta./(s + lambda(i)));
    eta(i) = norm(x_lambda(:,i));
  end
  if (nargout > 1 && size(U,1) > p)
    rho = sqrt(rho.^2 + norm(b - U(:,1:n)*[beta;U(:,p+1:n)'*b])^2);
  end

elseif (m>=n)

  % The overdetermined or square general-form case.
  x0 = V(:,p+1:n)*U(:,p+1:n)'*b;
  for i=1:ll
    xi = beta./(s(:,1) + lambda(i)*s(:,2));
    x_lambda(:,i) = V(:,1:p)*xi + x0;
    rho(i) = lambda(i)*norm(beta./(s(:,1)./s(:,2) + lambda(i)));
    eta(i) = norm(s(:,2).*xi);
  end
  if (nargout > 1 && size(U,1) > p)
    rho = sqrt(rho.^2 + norm(b - U(:,1:n)*[beta;U(:,p+1:n)'*b])^2);
  end

else

  % The underdetermined general-form case.
  x0 = V(:,p+1:m)*U(:,p+1:m)'*b;
  for i=1:ll
    xi = beta./(s(:,1) + lambda(i)*s(:,2));
    x_lambda(:,i) = V(:,1:p)*xi + x0;
    rho(i) = lambda(i)*norm(beta./(s(:,1)./s(:,2) + lambda(i)));
    eta(i) = norm(s(:,2).*xi);
  end

end
end