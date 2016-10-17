function [y_noise, Cn] = filtering (y, SNR, fs, model, f_low, f_high)
% function which adds noise, applies pre-filter and narrow-band filter to torso potentials. 
%
% This routine by Víctor Suárez Gutiérrez (victor.suarez.gutierrez@urjc.es)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUT:
% y: torso potentials.
% SNR: Signal to Noise Ratio
% fs: sampling frequency.
% model: place where rotor occurs {'SR', 'SAF', 'CAF'}
% f_low: low cutoff frequency for filtering.
% f_high: high cutoff frequency for filtering.
%
%OUTPUT:
% y_noise: filtered torso potentials.
% Cn: noise covariance matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng('default') % for reproducibility 

% % Add noise to y
[y_noise, Cn] = addwhitenoise(y, SNR); 

% Filter
y_noise = preprocessFA (y_noise, fs, f_low, f_high, model);

end   
        
function [ outSignal, Cn] = addwhitenoise( inSignal, SNR)
% Add gaussian white noise.

% We assume constant noise power in all electrodes, then different SNR for
% each electrodes. We consider 'SNR' the mean SNR value (in dB). 

% inSignal dimensions should be [M electrodes x T times];

PowerInSig   = mean( mean(abs(inSignal).^2,2) );   
PowerInSigdB = 10*log10(PowerInSig);  % dB units.

sigma     = sqrt(10^((PowerInSigdB-SNR) / 10 ));
noise = sigma*randn(size(inSignal));
outSignal = inSignal + noise;
Cn = (sigma^2)*eye(length(outSignal(:,1))); % noise covariance matrix.

end

