% Loads simuation parameters
% It also loads data set depending on simulation options.
%
% This routine by Víctor Suárez Gutiérrez (victor.suarez.gutierrez@urjc.es)
%
%
%% Paramters defined in this script:
% fs: sampling frecuency.
% f_low: low cut-off frecuency (default= 11 Hz).
% f_high: high cut-off frecuency (default= 13 Hz).
% params: paramters for each algorithm
% Mtransfer: transfer matrix.
% x: epicardial potentials [NxT]
% lambda_search_values: search grid for reg. param
% Nsamples: number of time instants


%% Common variables:
fs          = 500;       % Sampling frecuency.
Nsamples    = 2500;      % Number of instants.
Results     = [];        % Initialize.
params      = [];        % Initialize.
if strcmp(algorithm,'Tikhonov') || strcmp(algorithm,'TV')
    params = 10.^(-0.5:-0.5:-12.5); % Search grid for reg. param
end
SNR_name    = ['SNR' num2str(SNR)];

%% Low and high cutoff frequencies for bandpass filtering

if strcmp(model, 'SR')
    f_low = 0;
else
    f_low = 3;
end
f_high = 30;


%% Load algorithm's parameters
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
save_path = load_path;
    
if ~compute_params
    try
        params = load(params_file_name);
        params = params.param_opt;
    catch
        if ~compute_params,
            error('It is not possible to load parameters: file not found.');
        end

    end
end

%% Load geometry, transfer matrix and signals:
load([data_path 'geometry.mat']); 
load([data_path 'transfer.mat']);
load([data_path model '/' 'signal.mat']);

%% Remove samples during initial and end transients

switch model
    case 'SR'
        if strcmp(algorithm, 'Bayes')      % Leave 500 out for computing covariance matrix and mean
            Nsamples = 2000;        % Number of instants only for BayesMAP.
        end
        samples = 2:Nsamples;   %Sinusal sampling. Instant 1 is erased.
    case 'SAF'
        samples = length(EG(1,:))-Nsamples-1500+2:length(EG(1,:))-1500;
    case 'CAF'
        samples = length(EG(1,:))-Nsamples-2000+2:length(EG(1,:))-2000;
end
x = EG(:,samples);    % [NxT]; N: numbers of nodes (triangles´ vertexs) in epicardium.

% Window selection for Bayes method.
if strcmp(algorithm,'Bayes')
    if strcmp(model,'SR')
        ss = EG;
        ss(:,samples) = [];
        ss_random = datasample(ss(:,2:end),150,2,'Replace',false);
    end
    if strcmp(model,'SAF')
        ss = EG;
        ss = ss(:,500:999);
        ss_random = datasample(ss,150,2,'Replace',false);
    end
    if strcmp(model,'CAF')
        ss = EG;
        ss = ss(:,3500:4000);
        ss_random = datasample(ss,150,2,'Replace',false);
    end
end
clear EG;