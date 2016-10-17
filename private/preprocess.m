function [x,y,Cn,A] = preprocess (x, MTransfer, model, SNR, fs, f_low, f_high)
% Filtering epicarcial potentials, applying WCT correction and compute forward problem.
%
% This function by Víctor Suárez Gutiérrez (victor.suarez.gutierrez@urjc.es)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUT: 
% x: epicardial potentials [NxT]
% Mtransfer: transfer matrix.
% model: place where rotor occurs {'SR', 'SAF', 'CAF'}
% SNR: Signal to Noise Ratio.
% fs: sampling frecuency.
% f_low: low cut-off frecuency (default= 11 Hz).
% f_high: high cut-off frecuency (default= 13 Hz).
%
%OUTPUT:
% x: filtered epicardial potentials.
% y: filtered and noisy (SNR) torso potentials.
% Cn: noise covariance matrix.
% A: transfer matrix with WCT correction.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Remove mean and filter potentials
x = bsxfun(@minus,x,mean(x,2));  % Remove mean
x = preprocessFA (x, fs, f_low, f_high, model);

%% Transform transfer matrix to account for WCT correction
A_original = MTransfer;   % [MxN]; transfer matrix.
[n_torso,~] = size(A_original);
M_wct = (1/n_torso)*ones(n_torso);
A = (A_original - M_wct*A_original);    %A_wct

%% Forward problem
y = A*x;    % [MxT]; M: numbers of nodes (triangles´ vertexs) in torso.
[y, Cn] = filtering (y, SNR, fs, model, f_low, f_high);

