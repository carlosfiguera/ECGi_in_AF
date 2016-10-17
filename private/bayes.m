function x_hat = bayes (A, y, x, Cn, Nsamples, ss_random)
% % Bayesian MAP estimation method.
% 
% % This routine, by Carlos Figueras and Victor Suárez.
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Calculation:
xmean = mean(ss_random,2);
xmean = repmat(xmean, 1, length(x(1,:)));
Cx = (1/Nsamples) * (x-xmean)*(x-xmean).'; % Samples estimate of the covariance matrix
for j=1:length(y(1,:))
    fprintf('Time instant: %i\n',j);
    x_hat(:,j) = (Cx * A.' / (A*Cx*A.' + Cn)) * y(:,j);
end



