function [UI, OI, CC, MD] = drivermetrics (real, estimated, areas, xSMF, SMF, atrial_model)
% driver metrics calculation.
%
%by Víctor Suárez Gutiérrez (victor.suarez.gutierrez@urjc.es)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
% real: real modes.
% estimated: estimated modes.
% areas: area associated with nodes.
% xSMF: pdf does with real driver.
% SMF: pdf does with estimated driver.
% atrial_model: atria model.
%
% OUTPUT:
% UI: underestimation index.
% OI: overestimation index.
% CC: driver correlation coefficient.
% MD: euclidean distance between real and estimated drivers.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Assessment:
    nreal = numel(real(:,1));
    nestimated = numel(estimated(:,1));
    realTP = 0;
    estimatedTP = 0;
    FN = 0;
    FP = 0;
    for i=1:nreal
        if ~isempty(find(real(i,1)==estimated(:,1),1))
            realTP = realTP+real(i,2)*areas(real(i,1));
        else
            FN = FN+real(i,2)*areas(real(i,1));
        end
    end
    for i=1:nestimated
        if ~isempty(find(estimated(i,1)==real(:,1),1))
            estimatedTP = estimatedTP+estimated(i,2)*areas(estimated(i,1));
        else
            FP = FP+estimated(i,2)*areas(estimated(i,1));
        end
    end
                   
    % Driver metrics:
    UI = FN/(FN+realTP)*100; 
    OI = FP/(FP+estimatedTP)*100;
    realVALOR = zeros(size(areas(:,1)));
    realVALOR(real(:,1)) = real(:,2);
    VALOR = zeros(size(areas(:,1)));
    VALOR(estimated(:,1)) = estimated(:,2);
    [~, CC, ~] = timemetrics (realVALOR, VALOR);   
    
   % Euclidean distance between real and estimated drivers.
   MD = compute_SMF_distance(xSMF, SMF, atrial_model.distances);   
end


function [RMSE,CC,RDMS] = timemetrics (x,x_hat)
% This function calculates different metrics to evaluate the difference
% between real values (x) and estimated ones (x_hat)
%
% Input parameters:
%   - x: Ground Truth. Expected Dimensions N x T, 
%        being N the number of leads and T the number of time instants.
%   - x_hat: Estimatd signal. Expected dimensions N x T, 
%            being N the number of leads and T the number of time instants.
%
% This function by Felipe Alonso-Atienza (felipe.alonso@urjc.es)


% parsing
if nargin < 3
    verbose=0;
end

%% Norms
norm_x     = sqrt(sum(x.^2));
norm_x_hat = sqrt(sum(x_hat.^2));

%% Root mean square error
error = (x - x_hat);    % this is a matrix
RMSE = sqrt(sum(error.^2))./norm_x; % sum over rows.

%% Correlation coefficient
CC=sum(x.*x_hat)./(norm_x.*norm_x_hat);

%% Relative Difference Measure Star (RDMS)
% Definition taken from: INVERSE ELECTROCARDIOGRAPHY USING REDUCED 
% LEAD- SET BY TTLS AND LTTLS REGULARIZATION ALGORITHMS
%(http://www.wseas.us/e-library/conferences/2014/Istanbul/ELECT/ELECT-19.pdf)

% normalize signals by its power for each time instant
x_norm     = bsxfun(@rdivide,x,norm_x);         
x_norm_hat = bsxfun(@rdivide,x_hat,norm_x_hat);

% calculate the norm of the difference of signals (for each time instant)
difference = x_norm - x_norm_hat;
RDMS       = sqrt(sum(difference.^2));

%% mean values over time
MRMSE = mean(RMSE);
MCC   = mean(CC);

%% drawing for debugging
if verbose
    n = 0:numel(RMSE)-1;
    
    % RMSE and RDMS
    figure;
    plot(n, RMSE); 
    hold on;
    plot(n, RDMS, 'r'); 
    plot(n, ones(1, numel(n))*mean(RMSE), 'Linewidth', 4);
    plot(n, ones(1,numel(n))*mean(RDMS), 'r', 'Linewidth', 4);
    hold off;
    title('SNR 10 dB', 'Fontsize', 25);
    xlabel('n [sample]', 'Fontsize', 18);   
    ylabel('Normalized Error','Interpreter','Latex','Fontsize',18)
    legend('RMSE','RDMS')

    % CC
    figure;
    plot(n,CC); 
    hold on;
    plot(n,ones(1,numel(n))*mean(CC), 'Linewidth', 4);
    hold off;
    title('SNR 10 dB', 'Fontsize', 25);
    xlabel('n [sample]', 'Fontsize', 18);   
    ylabel('CC','Interpreter','Latex','Fontsize',18)
end

% MRMSE = mean(RMSE);
% MRDMS = mean(RDMS);
% MCC   = mean(CC);
end