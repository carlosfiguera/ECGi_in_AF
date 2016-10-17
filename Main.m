function Main()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         
% This code implements the experiments for the paper:
%   "Regularization Techniques for ECG Imaging during Atrial Fibrillation:
%   a Computational Study" 
%                                         
% Submitted to Frontiers in Physiology.
%
%               Felipe Alonso-Atienza & Víctor Suárez Gutiérrez &    
%               & Carlos Figuera Pozuelo &           
%
%               Sept. 2016
%                                                                         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all; close all; clc;

%% Options

model            = 'SAF'; % {'SR', 'SAF', 'CAF'}
SNR              = 20;    % (dBs)
compute_params   = 0;     % {1, 0} If 1, computes de regularization
                          % parameter and save to parameters.mat file. If
                          % 0, it loads de regularization parameter from
                          % parameters.mat file (stored in models folders).

algorithm        = 'Tikhonov';    % {'Tikhonov', 'TSVD', 'GS', 'TV',
                                 % 'Bayes', 'DSVD', 'GMRES'} -- Defined in
                                 % the paper.
                      
order            = 0;    % For Tikhonov and TSVD. Use [] for other methods.

reg_param_method = 'i';  % {'i','g'} Method for computing regularization 
                         % parameter: 'i': recompute for each instant;
                         % 'g': constant parameter (unique for all time
                         % instants. Use '' (empty string) for other
                         % methods.


%% Paths
data_path    = [pwd '/data/'];                       
results_path = [pwd '/results/'];
parameters_path = [pwd '/parameters/'];

%% Load main simulation parameters
load_parameters;

%% Pre-process potentials and compute forward problem
[x, y, Cn, A] = preprocess (x, MTransfer, model, SNR, fs, f_low, f_high);

%% Pre-compute DF, phase map and dirver for true potentials and A'A, L and L'L 
xDF = dominant_frequency (x, fs, model);   
xphase = instantphase (x, xDF, model, fs);
if ~strcmp(model, 'SR')
    [~, xmodes, xSMF] = driver_location (load([data_path 'filled_geometry.mat']) , x, xDF, fs);
end
[AA, L, LL] = precompute_matrices (A, order, atrial_model);

%% Epicardial potentials
param_opt = [];
switch algorithm
    
    case 'Tikhonov'
        [x_hat, param_opt] = tikhonov (A, L, AA, LL, y, params, SNR, order, reg_param_method, compute_params);
        
    case 'TSVD'
        [x_hat, param_opt] = tsvd (A, L, y, params, SNR, order, compute_params);
        
    case 'GS'
        [x_hat, param_opt] = greensite (A, L, y, params, SNR, compute_params);
        
    case 'TV'
        [x_hat, param_opt] = totalvariation (A, L, AA, LL, y, params);
        
    case 'Bayes'
        x_hat = bayes (A, y, x, Cn, Nsamples, ss_random);
        
    case 'DSVD'
        [x_hat, param_opt] = dsvd (A, y, params, SNR, compute_params);
          
    case 'GMRES'
        x_hat = gmresidual (A, y);
        
    otherwise
        error('Unknown algorithm');
        
end

% Dominant frequency
DF = dominant_frequency (x_hat, fs, model);
    
% Phase maps
phase = instantphase (x_hat, DF, model, fs);

% Driver Location
if ~strcmp(model, 'SR')
    [driver, modes, SMF] = driver_location (load([data_path 'filled_geometry.mat']) , x_hat, DF, fs);
end

%% Metrics
[Results] = metrics (Results, atrial_model, x, xDF, xphase, xmodes, xSMF, x_hat, DF, phase, driver, modes, SMF, SNR_name, [algorithm, num2str(order)], model, reg_param_method);

%% Save results
results_file = [results_path save_path];
save(results_file, 'Results', 'params', 'model', 'SNR', 'algorithm', 'order', 'reg_param_method', 'atrial_model', 'torso_model', 'x', 'y');

%% Save parameters
if compute_params
    parameters_file = [parameters_path save_path];
    save(parameters_file, 'param_opt');
end



