function visualize_restuls()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         
% This code shows the results of experiments for the paper:
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
SNR              = 20;   % (dBs)
algorithm        = 'Tikhonov';    % {'Tikhonov', 'TSVD', 'GS', 'TV',
                                 % 'Bayes', 'DSVD', 'GMRES'} -- Defined in
                                 % the paper.
                      
order            = 0;    % For Tikhonov and TSVD. Use [] for other methods.

reg_param_method = 'g';  % {'i','g'} Method for computing regularization 
                         % parameter: 'i': recompute for each instant;
                         % 'g': constant parameter (unique for all time
                         % instants. Use '' (empty string) for other
                         % methods.

flag             = 0;    % flag = 0 : torso representation (front and back sides).
                         % flag = 1 : atria potential representation of estimated potentials (both sides).
                         % flag = 2 : dominant frequencies representation of atria (both sides).    
                         % flag = 3 : representation of phases (both sides).
                         % flag = 4 : representation of drivers location and its pdf over the atria surface (both sides).
                         % flag = 5 : time metrics (both sides): RDMS and CC.
                         % flag = 6 : frequency metrics (both sides).
                         % flag = 7 : phase metrics (both sides): RDMS and CC.

video            = 0;     % get all frames when video=1 and only one instant when video is 0.


%% Paths
data_path    = [pwd '/data/'];                       
results_path = [pwd '/results/'];
parameters_path = [pwd '/parameters/'];
visualization_path = [pwd '/visualization/'];

%% Load algorithm's parameters
fs          = 500;       % Sampling frecuency.
Nsamples    = 2500;      % Number of instants.
Results     = [];        % Initialize.
params      = [];        % Initialize.
SNR_name    = ['SNR' num2str(SNR)];
switch algorithm
    case 'Tikhonov'
        load_path = ['/' model '/' algorithm '_' num2str(order) '_' reg_param_method '_' SNR_name];
        
    case 'TSVD'
        load_path = ['/' model '/' algorithm '_' num2str(order) '_' SNR_name];
        
    case 'GS'
        load_path = ['/' model '/' algorithm '_' SNR_name];
        
    case 'TV'
        load_path = ['/' model '/' algorithm '_' SNR_name];
        
    case 'Bayes'
        load_path = ['/' model '/' algorithm '_' SNR_name];
        
    case 'DSVD'
        load_path = ['/' model '/' algorithm '_' SNR_name];
        
    case 'GMRES'
        load_path = ['/' model '/' algorithm '_' SNR_name];
        
    otherwise
        error('Unknown algorithm');
end
params_file_name = [parameters_path load_path];
results_file_name = [results_path load_path];
save_path = load_path;
visualization_file_name = [visualization_path save_path];  

%% load Results and Parameters:
load(params_file_name);
load(results_file_name);

if isempty(Results)
    fprintf('There is no Results in these options.\n');
    clear all;
    return
end


%% Visualization of metrics:
visualization (atrial_model, torso_model, x, y, SNR, flag, video, model, algorithm, order, reg_param_method, Results, visualization_file_name);
