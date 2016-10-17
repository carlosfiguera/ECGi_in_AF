function phase = instantphase (ECG, DF, model, fs)
% Assessment of phase for node and instant.
%
% This function by Víctor Suárez Gutiérrez (victor.suarez.gutierrez@urjc.es)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUT: 
% ECG: signal.
% DF: dominant frequency.
% model: place where rotor occurs {'SR', 'SAF', 'CAF'}
% fs: sampling frequency.
%
%OUTPUT:
% phase: phases for each node.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nL,~]=size(ECG);
 
 if strcmp(model,'s')
     for i=1:nL
        ecg = preprocessFA (ECG(i,:), fs, 0, DF(i)+2, model);
        s=ecg-mean(ecg);
        ha=imag(hilbert(s));
        phase(i,:)=-atan2(ha,s);        
     end
 else
     for i=1:nL
        ecg = preprocessFA (ECG(i,:), fs, 3, DF(i)+2, model);
        s=ecg-mean(ecg);
        ha=imag(hilbert(s));
        phase(i,:)=-atan2(ha,s);        
     end
 end