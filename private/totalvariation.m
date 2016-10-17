function [x_hat, lambda_opt] = totalvariation (A, L, AA, LL, y, params)
% 
% % This routine, by Victor Suárez Gutiérrez.
% victor.suarez.gutierrez@urjc.es
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUT:
% A: transfer matrix.
% y: torso potentials.
% x: epicaldial potentials.
% Cn: noise covariance matrix.
% Nsamples: number of instances.
% ss_random: 150 random instants of epicardial potentials from a 1 second 
% window outside the window under experiment.

%OUTPUT:
% x_hat: estimated epicardial potentials.
% lambda_opt: optimal regularization parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for t=1:length(y)
    fprintf('Iteration number: %i.\n', t);
    y_t = y(:,t);
    for w=1:length(params)
        [~,cost,~] = lasso_admm_lap(y_t, A, L, AA, LL, params(w));   
        cost_params(w) = cost(end);
    end
    [~, bb] = min(cost_params);
    lambda_opt(1,t) = params(bb);
    [x_solution,~,~] = lasso_admm_lap(y_t, A, L, AA, LL, lambda_opt(1,t));
    x_hat(:,t) = x_solution;
end