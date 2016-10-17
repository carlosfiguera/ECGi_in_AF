function [x_hat, lambda_opt] = tikhonov (A, L, AA, LL, y, lambda, SNR, order, reg_param_method, compute_params)
% Tikhonov reconstruction method.
%
% This routine by Víctor Suárez Gutiérrez (victor.suarez.gutierrez@urjc.es)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUT:
% A: transfer matrix.
% L: matrix which takes part in regularization term of Tikhonov approach.
% AA: A'*A
% LL: L'*L
% y: torso potentials.
% lambda: regularization parameter.
% order: Tikhonov and TSVD order (0, 1 or 2).
% SNR: signal noise rate (dB).
% reg_param_method: global (g) or by instant (i) calculation.
% compute_params: 1 if reg_params need to be calculated. Otherwise 0. 

%
%OUTPUT:
% x_hat: estimated epicardial potentials.
% lambda_opt: optimal regularization parameters (one if global or an array by instant).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Determinate which method to be applied, by instants or global.
if (isempty(reg_param_method) || ~strcmp(reg_param_method,'g'))  % by instants.
    [x_hat, lambda_opt] = tikhonovinstants (y, A, L, lambda, SNR, order, compute_params);
else    %global
    [x_hat, lambda_opt] = Tikhonovglobalmethod (A, L, AA, LL, y, lambda, compute_params); 
end

end


function [x_hat, lambda_opt] = tikhonovinstants (y, A, L, lambda, SNR, order, compute_params)
%Tikhonov by instants reconstruction method.

% This routine by Víctor Suárez Gutiérrez (victor.suarez.gutierrez@urjc.es)

%%
if order~=0                            
    [U,s,V] = csvd(A,'full'); 
    [~, c, ~] = csvd(L,'full'); 
    c = c(1:length(s));
    s1 = s;
    s = s./c;
else
    [U,s,V] = csvd(A,'full');
    s1 = s;
end

if ~compute_params
    lambda_opt = lambda;
    for t=1:length(y(1,:))
        fprintf('Tikhonov Instant: t = %i.\n', t);
        x_hat(:,t) = V(:,1:numel(s))*diag((s.^2)./((s.^2)+lambda_opt(1,t)))*((U'*y(:,t))./s1);
    end
else
    for t=1:length(y(1,:))
        fprintf('Tikhonov Instant: t = %i.\n', t);
        [reg_corner, rho, eta, reg_param] = l_curve (U, s, s1, y(:,t), 'Tikh', SNR);
        lambda_opt(1,t) = reg_corner^2;
        x_hat(:,t) = V(:,1:numel(s))*diag((s.^2)./((s.^2)+lambda_opt(1,t)))*((U'*y(:,t))./s1);
    end    
end
end


function [x_hat, lambda_opt] = Tikhonovglobalmethod (A, L, AA, LL, y, lambda, compute_params)
%Tikhonov global method reconstruction.
%
% This routine by Víctor Suárez Gutiérrez (victor.suarez.gutierrez@urjc.es)
% using L-curve approach to Tikhonov Regularization. IEEE trans on biomed eng. 
% 2000. 47(9):1293-1296.
%
% The analytical solution of the inverse problem in terms of Tikhonov 
% regularization is 
% 
% phi_E_hat = inv(A'*A + lambda*L'*L)*A'*phi_T
%
%
% To reduce computational requirements we define
% AA = A'*A and LL = L'*L;

if ~compute_params
    invA = (AA+lambda*LL)\A.';
    x_hat = invA*y;
    lambda_opt = lambda;
else
    % we look for the opt value of Lambda
    % using the l-curve
        
    magnitude_term = zeros(size(lambda));
    error_term =  zeros(size(lambda));    
    
    for j=1:length(lambda)
        invA = (AA+lambda(j)*LL)\A.';
        x_hat = invA*y;
        error_term(j)    = norm(A*x_hat-y,'fro')^2;
        magnitude_term(j)= norm(L*x_hat,'fro')^2;
    end

    x = log10(error_term);
    z = log10(magnitude_term);
    dx = gradient(x,2); % uses h (x,h) as the spacing between points in each direction.
    dz = gradient(z,2); 
 
    % optimal lambda vector:
    ddx = gradient(dx,2);
    ddz = gradient(dz,2);
    curva =(dx.*ddz-ddx.*dz)./((dx.^2+dz.^2).^(3/2));
    maxcurva  = find(abs(curva) == max(abs(curva)), 1, 'first');
    
    lambda_opt = lambda(maxcurva);
    x_hat = ((AA+lambda_opt*LL)\A.')*y;   
    
    if 0     % change 0 for 1 for L-curve representation.
        figure;
        plot(x,z,'ko-');
        xlabel('log ||Ax-y||^2');
        ylabel('log ||Lx||^2');
        hold on;
        plot(x(maxcurva), z(maxcurva),'*r');
        hold on;
        plot(x(maxcurva+2), z(maxcurva+2),'*g');
        hold on;
        plot(x(maxcurva+1), z(maxcurva+1),'*');
    end
    
end
end
