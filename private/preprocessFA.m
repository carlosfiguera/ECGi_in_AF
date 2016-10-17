function [ECGP] = preprocessFA (ECG, fs, f_low, f_high, model)
% Frequency filtering of ECG between fpa-fpb.
% This routine by Víctor Suárez Gutiérrez (victor.suarez.gutierrez@urjc.es)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUT: 
% ECG: epicardial potentials [NxT]
% fs: sampling frecuency.
% f_low: low cut-off frecuency (default= 11 Hz).
% f_high: high cut-off frecuency (default= 13 Hz).
% model: place where rotor occurs {'SR', 'SAF', 'CAF'}
%
%OUTPUT:
% ECGP: filtered ECG.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Assessment:
[nL,~] = size(ECG);

if strcmp(model,'SR')
    for i=1:nL
        % Remove DC component
        ecg=ECG(i,:);
        ecg=ecg-mean(ecg);
        % LPF
        [b2,a2]=butter(4,f_high/round((fs/2)));   
        filteredsignal=filtfilt(b2,a2,ecg);
        ECGP(i,:)=filteredsignal;
    end
else
    for i=1:nL
        % Remove DC component
        ecg=ECG(i,:);
        ecg=ecg-mean(ecg);
        % HPF
        [b1,a1]=butter(4,f_low/(fs/2),'high');   
        filteredsignal=filtfilt(b1,a1,ecg);
        % LPF
        [b2,a2]=butter(4,f_high/(fs/2));   
        filteredsignal2=filtfilt(b2,a2,filteredsignal);
        ECGP(i,:)=filteredsignal2;
    end

end


        
