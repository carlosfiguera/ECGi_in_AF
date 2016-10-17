function [x_hat, k] = tsvd (A, L, y, k, SNR, order, compute_params)
% TSVD reconstruction method.
%
% This routine, by victor.suarez.gutierrez@urjc.es
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUT:
% A: transfer matrix.
% L: matrix which takes part in regularization term of Tikhonov approach.
% y: torso potentials.
% k: regularization parameter.
% SNR: signal noise rate (dB).
% order: Tikhonov order (0, 1 or 2).
% compute_params: True whether reg_params need to be calculated. Otherwise False. 

%OUTPUT:
% x_hat: estimated epicardial potentials.
% k: optimal regularization parameters (one if global or an array by instant).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if order~=0
        [U,sm,V] = cgsvd(A,L); %other option: [UU,sm,X] = cgsvd(A,L1);
    else
        [U,sm,V]   = csvd(A);
    end
        
    if ~compute_params
        for t=1:length(y(1,:))
            fprintf('TSVD Instant: t = %i.\n', t);
            if order==0
                x_hat(:,t) = tsvd0 (U, sm, V, y(:,t), k(:,t));
            else
                x_hat(:,t) = tgsvd (U, sm, V, y(:,t), k(:,t));
            end
        end
    else
        for t=1:length(y(1,:))
            fprintf('TSVD Instant: t = %i.\n', t);
            [reg_corner, rho, eta, reg_param] = l_curve (U, sm, sm, y(:,t), 'tsvd', SNR);
            k(:,t) = reg_corner;
            if order==0
                x_hat(:,t) = tsvd0 (U, sm, V, y(:,t), k(:,t));
            else
                x_hat(:,t) = tgsvd (U, sm, V, y(:,t), k(:,t));
            end
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


function [x_k,rho,eta] = tsvd0 (U,s,V,b,k)
%TSVD: Truncated SVD regularization.
%
% [x_k,rho,eta] = tsvdmethod(U,s,V,b,k)
%
% Computes the truncated SVD solution
%    x_k = V(:,1:k)*inv(diag(s(1:k)))*U(:,1:k)'*b .
% If k is a vector, then x_k is a matrix such that
%    x_k = [ x_k(1), x_k(2), ... ] .
% U, s, and V must be computed by the csvd function.
%
% The solution and residual norms are returned in eta and rho.

% Per Christian Hansen, DTU Compute, 12/21/97.

% Initialization.
[n,p] = size(V); lk = length(k);
if (min(k)<0 || max(k)>p)
  error('Illegal truncation parameter k')
end
x_k = zeros(n,lk);
eta = zeros(lk,1); rho = zeros(lk,1);
beta = U(:,1:p)'*b;
xi = beta./s;

% Treat each k separately.
for j=1:lk
  i = k(j);
  if (i>0)
    x_k(:,j) = V(:,1:i)*xi(1:i);
    eta(j) = norm(xi(1:i));
    rho(j) = norm(beta(i+1:p));
  end
end

if (nargout > 1 && size(U,1) > p)
  rho = sqrt(rho.^2 + norm(b - U(:,1:p)*beta)^2);
end
end


function [x_k,rho,eta] = tgsvd(U,sm,X,b,k)
%TGSVD Truncated GSVD regularization.
%
% [x_k,rho,eta] = tgsvd(U,sm,X,b,k) ,  sm = [sigma,mu]
%
% Computes the truncated GSVD solution
%            [ 0              0                 0    ]
%    x_k = X*[ 0  inv(diag(sigma(p-k+1:p)))     0    ]*U'*b .
%            [ 0              0             eye(n-p) ]
% If k is a vector, then x_k is a matrix such that
%    x_k = [ x_k(1), x_k(2), ... ] .
% U, sm, and X must be computed by the cgsvd function.
%
% The solution seminorm and the residual norm are returned in eta and rho.

% Reference: P. C. Hansen, "Regularization, GSVD and truncated GSVD",
% BIT 29 (1989), 491-504.

% Per Christian Hansen, DTU Compute, Feb. 24, 2008.

% Initialization
m = size(U,1);
n = size(X,1);
p = size(sm,1);
lk = length(k);
if (min(k)<0 | max(k)>p)
  error('Illegal truncation parameter k')
end
x_k = zeros(n,lk);
eta = zeros(lk,1); rho = zeros(lk,1);
beta = U'*b;
xi = beta(1:p)./sm(:,1);
if (nargout==3), mxi = sm(:,2).*xi; end

if (m>=n)

  % The overdetermined or square case. Treat each k separately.
  if (p==n)
    x_0 = zeros(n,1);
  else
    x_0 = X(:,p+1:n)*(U(:,p+1:n)'*b);
  end
  for j=1:lk
    i = k(j); pi1 = p-i+1;
    if(i==0)
      x_k(:,j) = x_0;
    else
      x_k(:,j) = X(:,pi1:p)*xi(pi1:p) + x_0;
    end
    if (nargout>1), rho(j) = norm(beta(1:p-i)); end
    if (nargout==3), eta(j) = norm(mxi(pi1:p)); end
  end

  if (nargout > 1 & size(U,1) > n)
    rho = sqrt(rho.^2 + norm(b - U(:,1:n)*beta(1:n))^2);
  end

else

  % The underdetermined case. Treat each k separately.
  if (p==m)
    x_0 = zeros(n,1);
  else
    x_0 = X(:,p+1:m)*(U(:,p+1:m)'*b);
  end
  for j=1:lk
    i = k(j); pi1 = p-i+1;
    if(i==0)
      x_k(:,j) = x_0;
    else
      x_k(:,j) = X(:,pi1:p)*xi(pi1:p) + x_0;
    end
    if (nargout>1), rho(j) = norm(beta(1:p-i)); end
    if (nargout==3), eta(j) = norm(mxi(pi1:p)); end
  end

end
end

function [U,s,V] = csvd(A,tst)
%CSVD Compact singular value decomposition.
%
% s = csvd(A)
% [U,s,V] = csvd(A)
% [U,s,V] = csvd(A,'full')
%
% Computes the compact form of the SVD of A:
%    A = U*diag(s)*V',
% where
%    U  is  m-by-min(m,n)
%    s  is  min(m,n)-by-1
%    V  is  n-by-min(m,n).
%
% If a second argument is present, the full U and V are returned.

% Per Christian Hansen, IMM, 06/22/93.

if (nargin==1)
  if (nargout > 1)
    [m,n] = size(A);
    if (m >= n)
      [U,s,V] = svd(full(A),0); s = diag(s);
    else
      [V,s,U] = svd(full(A)',0); s = diag(s);
    end
  else
    U = svd(full(A));
  end
else
  if (nargout > 1)
    [U,s,V] = svd(full(A)); s = diag(s);
  else
    U = svd(full(A));
  end
end
end


