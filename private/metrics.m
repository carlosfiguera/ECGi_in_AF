function [Results] = metrics (Results, atrial_model, x, xDF, xphase, xmodes, xSMF, x_hat, DF, phase, driver, modes, SMF, SNR, method, model, reg_param_method)
% Assessment of time, frequency, phase and driver metrics.
%
% This routine by Víctor Suárez Gutiérrez (victor.suarez.gutierrez@urjc.es)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%INPUT:
% Results: struc which stores final results.
% x: epicardial potentials.
% xDF: DF for real data in epicardium.
% xphase: phase of real epicadial potentials.
% xmodes: mode of each node where real driver appears {node mode}.
% xSMF: pdf does with real driver.
% x_hat: estimated epicardial potentials.
% DF; estimated DF.
% phase: estimated phase.
% driver: driver by instant.
% modes: mode of each node where estimtaed driver appears {node mode}.
% SMF: pdf does with estimated driver.
% SNR: Signal to Noise Ratio
% method: string defined by: ['method mode', num2str(order)]
% model: place where rotor occurs (SAF or CAF).
% reg_param_method: {'i','g'} method for computing regularization 
%
%OUTPUT:
% Results: struct which stores all results.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Time:
    [~, CC, RDMS] = timemetrics (x(1:2039,:), x_hat(1:2039,:));
    Results.(SNR).(model).(method).(reg_param_method).x_estimation = x_hat;
    Results.(SNR).(model).(method).(reg_param_method).RDMSn = RDMS';
    Results.(SNR).(model).(method).(reg_param_method).CCn = CC';
    [~, CC, RDMS] = timemetrics (x(1:2039,:)', x_hat(1:2039,:)');
    Results.(SNR).(model).(method).(reg_param_method).RDMSt = RDMS';
    Results.(SNR).(model).(method).(reg_param_method).CCt = CC';
    
    % DF:
    Results.(SNR).(model).(method).(reg_param_method).DF = DF(1:2039);

    % Frequency metric:
    Results.(SNR).(model).(method).(reg_param_method).RAE = abs(((xDF(1:2039)-DF(1:2039))./xDF(1:2039))*100);

    % Phase:
    [~, CC, RDMS] = timemetrics (xphase(1:2039,:)', phase(1:2039,:)');
    Results.(SNR).(model).(method).(reg_param_method).phase = phase;
    Results.(SNR).(model).(method).(reg_param_method).RDMSt_phase = RDMS';
    Results.(SNR).(model).(method).(reg_param_method).CCt_phase = CC';
        
    % Driver:
    [UI, OI, CCdriver, MD] = drivermetrics (xmodes, modes, atrial_model.areas, xSMF, SMF, atrial_model);
    Results.(SNR).(model).(method).(reg_param_method).drivers = driver;
    Results.(SNR).(model).(method).(reg_param_method).UI = UI;
    Results.(SNR).(model).(method).(reg_param_method).OI = OI;
    Results.(SNR).(model).(method).(reg_param_method).CCdriver = CCdriver;
    Results.(SNR).(model).(method).(reg_param_method).SMF = SMF;
    Results.(SNR).(model).(method).(reg_param_method).MD = MD;
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