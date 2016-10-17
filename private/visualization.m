function visualization (atrial_model, torso_model, x, y, SNR, flag, video, model, algorithm, order, reg_param_method, Results, visualization_file_name)
% Metrics and graphical results from different reconstruction models.
%
% This routine by Víctor Suárez Gutiérrez (victor.suarez.gutierrez@urjc.es)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUT:
% atrial_model: model of atria for 3D representations.
% torso_model: model of torso for 3D representations.
% x: real epicardial potentials.
% y: torso potencials.
% SNR: signal noise rate used to (dBs).
% flag: metrics which will be visualizated.
% video: 0=only a picture. 1=video.
% model: {'SR', 'SAF', 'CAF'}
% algorithm: {'Tikhonov', 'TSVD', 'GS', 'TV', 'Bayes', 'DSVD', 'GMRES'}
% order: for Tikhonov and TSVD. Use [] for other methods.
% reg_param_method:   {'i','g'} Method for computing regularization 
                    % parameter: 'i': recompute for each instant;
                    % 'g': constant parameter.
% Results: metrics from epicardial reconstructions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Load visualization parameters:
minF = 3;   % low frequency.
maxF = 8; % % high frequency.
if strcmp(model,'SR');
    instant_time = 826;
    minF = 1;   % low frequency.
    maxF = 1.4; % % high frequency.
end
if strcmp(model,'SAF');
    instant_time = 721;
end
if strcmp(model,'CAF');
    instant_time = 711;
end

%% Calculation:
if (video == 1) && ((flag == 2) || (flag == 5) || (flag == 6) || (flag == 7))
    fprintf('There is no video in this flag option.\n');
    return
end
if flag == 0
    normal = 10;    % dataset normalization.
    fprintf('Tikhonov SNR_%i:\n',SNR);
    if video
        my_frames = plotting (torso_model, y./normal, [], -1, 1, flag, video, [], SNR); 
        movie2avi(my_frames, [visualization_file_name, '_', 'flag' num2str(flag), '.avi']);  % save video.
        close all;
    else
        plotting (torso_model, y(:,instant_time)./normal, [], -1, 1, flag, video, instant_time, SNR); 
    end
end

if flag == 1
    normal = 50;
    resultset = Results.(['SNR' num2str(SNR)]).(model).([algorithm num2str(order)]).(reg_param_method).x_estimation;
    resultset = resultset(1:2039,:)./normal;
    fprintf([algorithm ' SNR_%i:\n'],SNR);
    if video
        my_frames = plotting (atrial_model, resultset, [], -1, 1, flag, video, [], SNR); 
        movie2avi(my_frames, [visualization_file_name, '_', 'flag' num2str(flag), '.avi']);  % save video.        
        close all;
    else
        plotting (atrial_model, resultset(:,instant_time), [], -1, 1, flag, video, instant_time, SNR); 
    end
end

if flag == 2
    video = 0;
    resultset = Results.(['SNR' num2str(SNR)]).(model).([algorithm num2str(order)]).(reg_param_method).DF;
    resultset = resultset(1:2039,:);
    fprintf([algorithm ' SNR_%i:\n'],SNR);
    plotting (atrial_model, resultset(:,1), [], minF, maxF, flag, video, instant_time, SNR); 
end

if flag == 3
    resultset = Results.(['SNR' num2str(SNR)]).(model).([algorithm num2str(order)]).(reg_param_method).phase;
    resultset = resultset(1:2039,:);
    fprintf([algorithm ' SNR_%i:\n'],SNR);
    if video
        my_frames = plotting (atrial_model, resultset, [], -pi, pi, flag, video, [], SNR); 
        movie2avi(my_frames, [visualization_file_name, '_', 'flag' num2str(flag), '.avi']);  % save video.
        close all;
    else
        plotting (atrial_model, resultset(:,instant_time), [], -pi, pi, flag, video, instant_time, SNR); 
    end
end

if (flag == 4) && (~strcmp(model, 'SR'))
    fprintf([algorithm ' SNR_%i:\n'],SNR);
    if video
        drivers = Results.(['SNR' num2str(SNR)]).(model).([algorithm num2str(order)]).(reg_param_method).drivers;
        drivers2 = zeros(2039,length(drivers(:,1)));
        for i = 1:length(drivers2(1,:))
            if sum(drivers(i,:)) ~= 0
                drivers2(drivers(i,drivers(i,:)~=0),i) = 1;
            end
        end
        drivers = drivers2;
        clear drivers2;
        my_frames = plotting (atrial_model, drivers, [], 0, 1, flag, video, [], SNR); 
        movie2avi(my_frames, [visualization_file_name, '_', 'flag' num2str(flag), '.avi']);  % save video.
        close all;
    else
        SMF = Results.(['SNR' num2str(SNR)]).(model).([algorithm num2str(order)]).(reg_param_method).SMF;
        if strcmp(model, 'SAF')
            plotting (atrial_model, SMF, [], 0, 0.1, flag, video, instant_time, SNR); 
        else
            plotting (atrial_model, SMF, [], 0, 0.4, flag, video, instant_time, SNR); 
        end
    end
end

if flag == 5
    video = 0;
    RDMSt = Results.(['SNR' num2str(SNR)]).(model).([algorithm num2str(order)]).(reg_param_method).RDMSt;
    RDMSt = RDMSt(1:2039,:);
    CCt = Results.(['SNR' num2str(SNR)]).(model).([algorithm num2str(order)]).(reg_param_method).CCt;
    CCt = CCt(1:2039,:);
    fprintf([algorithm ' SNR_%i:\n'],SNR);
    plotting (atrial_model, RDMSt, [], min(RDMSt), max(RDMSt), flag, video, instant_time, SNR); 
    figure;
    plotting (atrial_model, CCt, [], min(CCt), max(CCt), flag, video, instant_time, SNR); 
end

if flag == 6
    video = 0;
    RAE = Results.(['SNR' num2str(SNR)]).(model).([algorithm num2str(order)]).(reg_param_method).RAE;
    fprintf([algorithm ' SNR_%i:\n'],SNR);
    plotting (atrial_model, RAE, [], 0.25*mean(RAE), 4*mean(RAE), flag, video, instant_time, SNR); 
end

if flag == 7
    video = 0;
    RDMSt_phase = Results.(['SNR' num2str(SNR)]).(model).([algorithm num2str(order)]).(reg_param_method).RDMSt_phase;
    RDMSt_phase = RDMSt_phase(1:2039,:);
    CCt_phase = Results.(['SNR' num2str(SNR)]).(model).([algorithm num2str(order)]).(reg_param_method).CCt_phase;
    CCt_phase = CCt_phase(1:2039,:);
    fprintf([algorithm ' SNR_%i:\n'],SNR);
    plotting (atrial_model, RDMSt_phase, [], min(RDMSt_phase), max(RDMSt_phase), flag, video, instant_time, SNR); 
    figure;
    plotting (atrial_model, CCt_phase, [], min(CCt_phase), max(CCt_phase), flag, video, instant_time, SNR); 
end

%% Save figure
if video == 0
    nfigs = length(get(0,'Children'));
    if nfigs == 1
        savefig([visualization_file_name 'flag' num2str(flag)])
        close all;
    else
        metric_name = {'CC', 'RDMS'};
        for i=1:nfigs
            get(groot,'CurrentFigure');
            savefig([visualization_file_name 'flag' num2str(flag) '_' metric_name{i}])
            close gcf;
        end
    end
end
end


function my_frames = plotting (Model, MAG, EST, ming, maxg, flag, video, instant_time, SNR)
% It displays a video or a picture.

% Input variables:
% Model: structure´s model.
% MAG: data of the structure´s model.
% EST: data of the estimated signal.
% ming: minimal value of the scale.
% maxg: maximum value of the scale.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if video
    [~, c] = size(MAG);
    for i=1:c
        if flag == 0
            cmap = buildcmap('bcyr');
            colormap(cmap);
            aux = torsorepresentation (Model, MAG(:,i), ming, maxg);
            
            ax = axes('Units','Normal','Position',[.075 .075 .87 .82],'Visible','off');
            set(get(ax,'Title'),'Visible','on')
            text = (['Time :  t = ' num2str(i)]);
            title(text, 'Fontsize',25);

            ax = axes('Units','Normal','Position',[.075 .075 .88 .73],'Visible','off');
            set(get(ax,'Title'),'Visible','on')
            text = (['SNR ' num2str(SNR)]);
            title(text, 'Fontsize',18);

            ax = axes('Units','Normal','Position',[.075 .075 .56 .68],'Visible','off');
            set(get(ax,'Title'),'Visible','on')
            text = ('Front');
            title(text, 'Fontsize',14);

            ax = axes('Units','Normal','Position',[.075 .075 1.19 .68],'Visible','off');
            set(get(ax,'Title'),'Visible','on')
            text = ('Back');
            title(text, 'Fontsize',14);
        end
        if (flag == 1) || (flag == 3) || (flag == 4)
            if flag == 1
                cmap = buildcmap('bcyr');
                colormap(cmap);
            end
            if flag == 3
                cmap = buildcmap('wyrgbmk');
                colormap(cmap);
            end
            if flag == 4
                cmap = buildcmap('bwr');
                cmap=cmap(30:end,:);
                colormap(cmap);
            end
            aux = atrialrepresentation (Model, MAG(:,i), ming, maxg);
            
            ax = axes('Units','Normal','Position',[.075 .075 .82 .82],'Visible','off');
            set(get(ax,'Title'),'Visible','on')
            text = (['Time :  t = ' num2str(i)]);
            title(text, 'Fontsize',25);

            ax = axes('Units','Normal','Position',[.075 .075 .83 .73],'Visible','off');
            set(get(ax,'Title'),'Visible','on')
            text = (['SNR ' num2str(SNR)]);
            title(text, 'Fontsize',18);

            ax = axes('Units','Normal','Position',[.075 .075 .46 .65],'Visible','off');
            set(get(ax,'Title'),'Visible','on')
            text = ('Front');
            title(text, 'Fontsize',14);

            ax = axes('Units','Normal','Position',[.075 .075 1.19 .65],'Visible','off');
            set(get(ax,'Title'),'Visible','on')
            text = ('Back');
            title(text, 'Fontsize',14);
            clear colormap;
        end
        pause(0.0005)

        my_frames(i) = getframe(gcf);   % frame saved.
        clf
    end
else
    my_frames = [];
    if flag == 0
        cmap = buildcmap('bcyr');
        colormap(cmap);
        aux = torsorepresentation (Model, MAG, ming, maxg);
        
        ax = axes('Units','Normal','Position',[.075 .075 .87 .85],'Visible','off');
        set(get(ax,'Title'),'Visible','on')
        text = (['Time :  t = ' num2str(instant_time)]);
        title(text, 'Fontsize',25);

        ax = axes('Units','Normal','Position',[.075 .075 .88 .80],'Visible','off');
        set(get(ax,'Title'),'Visible','on')
        text = (['SNR ' num2str(SNR)]);
        title(text, 'Fontsize',22);

        ax = axes('Units','Normal','Position',[.075 .075 .56 .79],'Visible','off');
        set(get(ax,'Title'),'Visible','on')
        text = ('Front');
        title(text, 'Fontsize',18);

        ax = axes('Units','Normal','Position',[.075 .075 1.19 .79],'Visible','off');
        set(get(ax,'Title'),'Visible','on')
        text = ('Back');
        title(text, 'Fontsize',18);
    end
    if (flag == 1) || (flag == 3)
        if flag == 1
            cmap = buildcmap('bcyr');
            colormap(cmap);
        end
        if flag == 3
            cmap = buildcmap('wyrgbmk');
            colormap(cmap);
        end
        aux = atrialrepresentation (Model, MAG, ming, maxg);
        ax = axes('Units','Normal','Position',[.075 .075 .81 .85],'Visible','off');
        set(get(ax,'Title'),'Visible','on')
        text = (['Time :  t = ' num2str(instant_time)]);
        title(text, 'Fontsize',25);

        ax = axes('Units','Normal','Position',[.075 .075 .815 .80],'Visible','off');
        set(get(ax,'Title'),'Visible','on')
        text = (['SNR' num2str(SNR)]);
        title(text, 'Fontsize',18);
    
        ax = axes('Units','Normal','Position',[.075 .075 .46 .75],'Visible','off');
        set(get(ax,'Title'),'Visible','on')
        text = ('Front');
        title(text, 'Fontsize',18);

        ax = axes('Units','Normal','Position',[.075 .075 1.19 .75],'Visible','off');
        set(get(ax,'Title'),'Visible','on')
        text = ('Back');
        title(text, 'Fontsize',18);
    end
    if (flag == 2) || (flag == 4) || (flag == 5) || (flag == 6) || (flag == 7)
        if flag == 4
            cmap = buildcmap('bwr');
            cmap=cmap(30:end,:);
            colormap(cmap);
        end
        aux = atrialrepresentation (Model, MAG, ming, maxg);
        ax = axes('Units','Normal','Position',[.075 .075 .815 .80],'Visible','off');
        set(get(ax,'Title'),'Visible','on')
        text = (['SNR' num2str(SNR)]);
        title(text, 'Fontsize',18);
    
        ax = axes('Units','Normal','Position',[.075 .075 .46 .75],'Visible','off');
        set(get(ax,'Title'),'Visible','on')
        text = ('Front');
        title(text, 'Fontsize',18);

        ax = axes('Units','Normal','Position',[.075 .075 1.19 .75],'Visible','off');
        set(get(ax,'Title'),'Visible','on')
        text = ('Back');
        title(text, 'Fontsize',18);
        clear colormap;
    end
end

end


function aux = torsorepresentation (Model, MAG, ming, maxg)
% It displays an instant time.

% Input variables:
% Model: structure´s torso.
% MAG: potentials.
% ming: minimal value of the scale.
% maxg: maximum value of the scale.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

faces    = Model.faces; % in this structure´s model they are triangles.
vertices = Model.vertices;  % nodes.

if nargin >2
    caxis([ming maxg]);
else
    maxCAXIS=max(max(abs( MAG(:,2:end) )));
    caxis(maxCAXIS.*[-1 1]);
end

window1=subplot(121);
aux = patch('Faces', faces,'Vertices', vertices,'FaceVertexCData', MAG); 
xlabel('X'); ylabel('Y'); zlabel('Z');
axis equal;

shading faceted;
axis off
caxis([ming maxg]);

% Lighting
lighting phong
shading interp;
view(0,0); 
zoom(2);
set(window1,'Position',[0.15 0.05 0.42 1]);

window2=subplot(122);
eaux = patch('Faces', faces,'Vertices', vertices,'FaceVertexCData', MAG);
axis equal;
shading faceted;
axis off
caxis([ming maxg]);

% Lighting
lighting phong
shading interp;
view(180,0); 
set(window2,'Position',[0.525 0 0.278 1.098]);
h=colorbar('southoutside');
set(h, 'Position', [0.2 0.2 .6 .03]);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',45)
end


function aux = atrialrepresentation (Model, MAG, ming, maxg)
% It displays an instant time.

% Input variables:
% Model: structure´s atria.
% MAG: potentials.
% ming: minimal value of the scale.
% maxg: maximum value of the scale.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

faces    = Model.faces; % in this structure´s model they are triangles.
vertices = Model.vertices;  % nodes.

if nargin >2
    caxis([ming maxg]);
else
    maxCAXIS=max(max(abs( MAG(:,2:end) )));
    caxis(maxCAXIS.*[-1 1]);
end

window1=subplot(121);
aux = patch('Faces', faces,'Vertices', vertices,'FaceVertexCData', MAG);  
xlabel('X'); ylabel('Y'); zlabel('Z');
axis equal;
shading faceted;
axis off
caxis([ming maxg]);

% Lighting
lighting phong
shading interp;
view(90,30);
set(window1,'Position',[0.09 0.05 0.42 1]);

window2=subplot(122);
eaux = patch('Faces', faces,'Vertices', vertices,'FaceVertexCData', MAG);   
xlabel('X'); ylabel('Y'); zlabel('Z');
axis equal;
shading faceted;
axis off
caxis([ming maxg]);

% Lighting
lighting phong
shading interp;
view(-90,30);
zoom(1.13);
set(window2,'Position',[0.48 0 0.42 1]);

h=colorbar('southoutside');
set(h, 'Position', [0.1 0.2 .815 .03]);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',45)

end


