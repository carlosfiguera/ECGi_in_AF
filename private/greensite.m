function [X_reg, lambda_opt] = greensite (A, L, y, lambda, SNR, compute_params)
% Greensite reconstruction method.
%
% This function by Víctor Suárez Gutiérrez (victor.suarez.gutierrez@urjc.es)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUT:
% A: transfer matrix.
% L: matrix which takes part in regularization term of Tikhonov approach.
% y: torso potentials.
% lambda: regularization parameter.
% SNR: signal noise rate (dB).
% compute_params: True whether reg_params need to be calculated. Otherwise False. 

%OUTPUT:
% X_reg: estimated epicardial potentials.
% lambda_opt: optimal regularization parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[U,~,~] = svd(y'*y);
trans_Y=y*U;
[Ua,sa,Va] = csvd(A,'full'); 
[~, c, ~] = csvd(L,'full'); 
c = c(1:length(sa));
s1 = sa;
sa = sa./c;

if ~compute_params
    lambda_opt = lambda;
    for t=1:length(trans_Y(1,:))
        fprintf('GS Instant: t = %i.\n', t);
        X_reg(:,t) = Va(:,1:numel(sa))*diag((sa.^2)./((sa.^2)+lambda_opt(1,t)))*((Ua'*trans_Y(:,t))./s1);
    end
else
    for t=1:length(trans_Y(1,:))
        fprintf('GS Instant: t = %i.\n', t);
        [reg_corner, rho, eta, reg_param] = l_curve (Ua, sa, s1, trans_Y(:,t), 'grst', SNR);
        lambda_opt(1,t) = reg_corner.^2;
        X_reg(:,t) = Va(:,1:numel(sa))*diag((sa.^2)./((sa.^2)+lambda_opt(1,t)))*((Ua'*trans_Y(:,t))./s1);
    end
end
end




  
  
 
    
