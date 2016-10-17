function xDF = dominant_frequency (x, fs, model)
% Assessment of dominant frequency for each node.
%
% This function by Víctor Suárez Gutiérrez (victor.suarez.gutierrez@urjc.es)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUT: 
% x: signal.
% fs: sampling frequency.
% model: place where rotor occurs {'SR', 'SAF', 'CAF'}
%
%OUTPUT:
% xDF: dominant frequency for each node of model.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Assessment:
[N,~] = size(x);
xDF = zeros(N,1);

if strcmp(model,'SR')
    for i=1:N
        signal = x(i,:);
        rx = xcorr(signal,'coeff');
        [pks,locs] = findpeaks(rx,'MinPeakDistance',200);
        xDF(i,:) = fs*(round(length(locs)/2)-1)/length(signal(1,:)); 
    end
else    
    L = 1024;
    j=1;
    while L > length(x(1,:))
        L = 2^nextpow2(N-j);
        j=j+1;
    end 
    for i=1:N
        egm  = x(i,:);
        egm = egm -mean(egm);    
        [Pff,f] = pwelch(egm,hamming(L),[],L,fs);
        [~,idx] = max(Pff);
        % Multiple of 3:
        [~,pf] = min(abs(f-f(idx)/3));
        thmin = f(pf)-0.5;
        thmax = f(pf)+0.5;
        [~,idxmin1] = min(abs(f-thmin));
        [~,idxmax1] = min(abs(f-thmax));
        % Multiple of 2:
        [~,pf] = min(abs(f-f(idx)/2));
        thmin = f(pf)-0.5;
        thmax = f(pf)+0.5;
        [~,idxmin2] = min(abs(f-thmin));
        [~,idxmax2] = min(abs(f-thmax));
        if ~isempty(find(Pff(idxmin1:idxmax1)>0.1*max(Pff))) 
            [~,pos] = max(Pff(idxmin1:idxmax1));
            pos = pos+idxmin1-1;
        elseif ~isempty(find(Pff(idxmin2:idxmax2)>0.1*max(Pff)))
            [~,pos] = max(Pff(idxmin2:idxmax2));
            pos = pos+idxmin2-1;
        else
            [~, pos] = max(Pff);
        end

        xDF(i) = f(pos);
    end
end