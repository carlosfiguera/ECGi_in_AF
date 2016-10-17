function [finalrotors, modes, SMF] = driver_location (filled_geometry, x, DF, fs)
% Assessment of driver location and its probability density function.
%
% This routine by Víctor Suárez Gutiérrez (victor.suarez.gutierrez@urjc.es)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUT:
% filled_geometry: model of filled atria for 3D representations.
% x: real epicardial potentials.
% DF: Dominant Frequency.
% fs: sample frequency.

%OUTPUT:
% finalrotors: driver location at each instant.
% modes: driver mode over nodes of atria surface.
% SMF: pdf of driver location.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

results = driver_searching (filled_geometry.atrial_model,x',fs);
drivers = results.pos_rot;
drivers(drivers>2039)=0;
freqdetec = unique(DF);
for i=1:length(freqdetec)
    moda(i,1)=freqdetec(i);
    moda(i,2) = sum(DF==moda(i,1));
end
f_rotor = max(moda(moda(:,2)>=100));

rotlocation = find(DF==f_rotor);
finalrotors = [];
ref = 1.2;   % cycles for being considerated as rotor.
threshold1 = round(fs*(ref/f_rotor));  % samples of rotor. time  = 1000*(ref/f_real_rotor);
for i=1:length(drivers(1,:))
    if (~isempty(intersect(drivers(:,i)',rotlocation')))
        if numel(find(drivers(:,i))) > threshold1
            finalrotors = [finalrotors drivers(:,i)];
        end
    end
end
drivers = finalrotors;
finalrotors = [];
for i=1:length(drivers(1,:))
    if (~isempty(intersect(drivers(:,i)',rotlocation')))
            drivers(find(~ismember(drivers(:,i),rotlocation)),i)=0;
            finalrotors = [finalrotors drivers(:,i)];
    end
end

%mode of each node:
modes = unique(finalrotors);
modes = modes(modes~=0);
for i=1:length(modes)
    modes(i,2) = sum(sum(finalrotors==modes(i,1)));
end
SMF = zeros(size(x(1:2039,1)));
pdf=modes(:,2)/sum(modes(:,2));
SMF(modes(:,1)) = pdf;

