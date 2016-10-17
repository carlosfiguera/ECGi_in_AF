function [x_filt, AA, L, LL] = preprocessing (x, fs, fc1, fc2, A, F, order, Aur)
% Function which applies narrow-band filter to x and calculates variables for 
% minimizing computing time in calculation of inverse matrix A.
%
% This routine by Víctor Suárez Gutiérrez (victor.suarez.gutierrez@urjc.es)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUT:
% x: epicardial potentials.
% fs: sampling frequency.
% fc1: low cutoff frequency for narrow-band filter.
% fc2: high cutoff frequency for narrow-band filter.
% A: transfer matrix.
% F: it is 1 if narrow-band filter is applied. Otherwise 0.
% order: Tikhonov order (0, 1 or 2).
% Aur: atrial model.
%
%OUTPUT:
% x_filt: filtered signal of epicardial potentials.
% AA: A'*A
% L: matrix which takes part in regularization term of Tikhonov approach.
% LL: L'*L
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% narrow-band filter:
    if (strcmp(F,'Yes') || strcmp(F,'Y') || strcmp(F,'y'))
        [x_filt, ~] = filtrarenbandaestrecha (x, fs, [fc1 fc2]);
    else
        x_filt = [];
    end

%% minimizing computing time:
    if ~isnan(order)
        [AA, L, LL] = variables_tik (A, order, Aur);
    else
        AA = [];
        L = [];
        LL = [];
    end
end


function [ECGf,fsd]=filtrarenbandaestrecha(ECG,fs,frecuencias)

% ECG tiene todas las derivaciones de ECG por fila
% fs es la frecuencia de muestreo
% frecuencias es un vector de dos columnas
%                    1ª columna: frecuencias de corte inferior
%                    2ª columna: frecuencias de corte superior
% ECGf contiene las derivaciones filtradas
% fsd es la frecuencia de muestreo de salida


pintarfiltros=0;
[nL,nT]=size(ECG);

% Calculo los parametros de diezmado para que sea un numero entero cercano a 128 Hz

ND=round(fs/128);    %coeficiente de diezmado
fsd=fs./ND;

% Primero diseño el paso bajo
d1=fdesign.lowpass;
setspecs(d1,'N,Fp,Ap',10,frecuencias(2),0.01,fsd);
hd = cheby1(d1);
[b,a]=sos2tf(hd.sosMatrix,hd.scaleValues);
fpb.b=b;
fpb.a=a;

% Y despues el paso alto
d2=fdesign.highpass;
setspecs(d2,'N,Fp,Ap',10,frecuencias(1),0.01,fsd);
hd = cheby1(d2);
[b,a]=sos2tf(hd.sosMatrix,hd.scaleValues);    
fpa.b=b;
fpa.a=a;
    
% Representacion de las respuestas en frecuencia de los filtros
if pintarfiltros==1
    [h1,w1] = freqz(fpb.b,fpb.a,fsd*2,fsd);
    [h2,w2] = freqz(fpa.b,fpa.a,fsd*2,fsd);
    %figure;
    plot(w1,20.*log10(h1));
    hold on; plot(w2,20.*log10(h2),'g');
end
    

% Para inicializar los datos de salida primero diezmo una derivacion
t = (0:nT-1)/fs;

der1_dec=decimate(ECG(1,:),ND);
nSd=length(der1_dec);
ECGf=zeros(nL,nSd);
%ECGfi=zeros(nL,nT);

%keyboard;

for i=1:nL
    sig_dec=decimate(ECG(i,:),ND);
    sig_fpb=filtfilt(fpb.b,fpb.a,sig_dec);
    ECGf(i,:)=filtfilt(fpa.b,fpa.a,sig_fpb);
    
    %ECGfi(i,:) = interpl(ECGf(i,:),t(1:ND:end),t);
    
end
end


function [AA, L, LL] = variables_tik (A, order, Modelo)
%Routine for calculating variables used to Tikhonov approaches: A'*A y L'*L.

% This routine by Víctor Suárez Gutiérrez (victor.suarez.gutierrez@urjc.es)

%%
switch order
    case 0,     
        [~,T] = size(A);
        L = eye(T); % identity matrix.
    case 1,
        xv  = Modelo.vertices(:,1);
        yv  = Modelo.vertices(:,2);
        zv  = Modelo.vertices(:,3);
        tri = Modelo.faces; 
        L = gradiente_mat(xv,yv,zv,tri,A);
    case 2,
        xv  = Modelo.vertices(:,1);
        yv  = Modelo.vertices(:,2);
        zv  = Modelo.vertices(:,3);
        tri = Modelo.faces; 
        L = laplaciano_mat(xv,yv,zv,tri,A);
end

AA = A'*A;
LL = L'*L; % norm 2

end


function [R,GRAD] = gradiente_mat(x,y,z,elementos,F)
% Assessment of gradient.

% APROXIMACIÓN 2 AL GRADIENTE
% esta funcion calcula el gradiente de una funcion en una superficie de dos
% dimensiones que vive en tres. Se supone que dicha superficie esta
% definida como una malla triangular de nodos. x y z son las componentes de
% cada nodo. En elementos estan puesto los nodos que forman cada uno de los
% elemento. Suponemos trianguares. No se ha tenido en cuenta si la
% superficie posee frontera !!!Cuidadín¡¡¡¡

% Igual que el laplaciano pero sin dividir por el área

%Inicializar variables:
[~,N] = size(F);
% GRAD  = zeros(size(F));
R     = zeros(N,N);

[num_ele aux]   = size(elementos);
num_nodos       = size(x);
area_elemento_i = zeros(num_ele,1);
vecinos=recopila(elementos);    %obtiene el número de faces al que pertenece un nodo y el índice de dichos faces.

% calcular el área de cada triángulo
for elemento=1:num_ele
    n_e = elementos(elemento,:);    %nodos ubicados en los vértices.
    P1=[x(n_e(1)),y(n_e(1)),z(n_e(1))];     % coordenadas de cada nodo.
    P2=[x(n_e(2)),y(n_e(2)),z(n_e(2))];
    P3=[x(n_e(3)),y(n_e(3)),z(n_e(3))];
    area_elemento_i(elemento)=area_triangulo(P1,P2,P3); %área del triángulo.
end
    
% Construir la matriz B
% para cada nodo 
for i_nodo=1:num_nodos
    area_total=0; 
    num_vecinos=vecinos(i_nodo,1);

    % por cada elemento que contenga al nodo (vecino)
    for i_vecino=1:num_vecinos;

        vecino_de_i=vecinos(i_nodo,i_vecino+1);
        i_nodo_vecino=elementos(vecino_de_i,:); %nodos que forman face 'vecino_de_i'.

        for nodos=1:3 

            % ii, jj, kk interesantes: índice de los distintos nodos
            % del elemento (ii el nodo en cuestión)
            ii=i_nodo;                  
            jj=i_nodo_vecino(nodos);     
            for nod=1:3
                if i_nodo_vecino(nod)~=ii && i_nodo_vecino(nod)~=jj
                    kk=i_nodo_vecino(nod);
                end
            end               

            % r_i, r_j, r_k: coordenadas de ii,jj,kk
            r_i=[x(ii),y(ii),z(ii)];
            r_j=[x(jj),y(jj),z(jj)]; 
            r_k=[x(kk),y(kk),z(kk)];  

            % es por la contribución doble que está el 1/2? -> por qué hace
            % esto?
            if ii~=jj
                R(i_nodo,jj) = R(i_nodo,jj)+(1/2)*(norm(r_j-r_k,2))/ norm(r_i-r_j,2);
                R(i_nodo,ii) = R(i_nodo,ii) - ((1/2)*(norm(r_j-r_k,2))/norm(r_i-r_j,2));
            end 
        end 

        % area sumada de todos los vecinos (elementos que contienen a un
        % nodo)
        area_total = area_total + area_elemento_i(vecino_de_i);
    end 

    %B(i_nodo,:) = (1/area_total)*B(i_nodo,:);  
end  

% GRAD = R*F;
GRAD = 0;
end


function [R,LAP]=laplaciano_mat(x,y,z,elementos,F)
% Assessment of laplacian.

% APROXIMACIÓN 2 AL LAPLACIANO
% esta funcion calcula el laplacino de una funcion en una superficie de dos
% dimensiones que vive en tres. Se supone que dicha superficie esta
% definida como una malla triangular de nodos. x y z son las componentes de
% cada nodo. En elementos estan puesto los nodos que forman cada uno de los
% elemento. Suponemos trianguares. No se ha tenido en cuenta si la
% superficie posee frontera !!!Cuidadín¡¡¡¡

[~,N] = size(F);
% LAP = zeros(size(F));
R = zeros(N,N);

[num_ele aux] = size(elementos);
num_nodos = size(x);
area_elemento_i = zeros(num_ele,1);
vecinos=recopila(elementos);

% calcular el área de cada triángulo
for elemento=1:num_ele
    n_e=elementos(elemento,:);
    P1=[x(n_e(1)),y(n_e(1)),z(n_e(1))];
    P2=[x(n_e(2)),y(n_e(2)),z(n_e(2))];
    P3=[x(n_e(3)),y(n_e(3)),z(n_e(3))];
    area_elemento_i(elemento)=area_triangulo(P1,P2,P3);
end
    
% Construir la matriz B
% para cada nodo 
for i_nodo=1:num_nodos
    area_total=0; 
    num_vecinos=vecinos(i_nodo,1);

    % por cada elemento que contenga al nodo (vecino)
    for i_vecino=1:num_vecinos;

        vecino_de_i=vecinos(i_nodo,i_vecino+1);
        i_nodo_vecino=elementos(vecino_de_i,:);

        for nodos=1:3 

            % ii, jj, kk interesantes: índice de los distintos nodos
            % del elemento (ii el nodo en cuestión)
            ii=i_nodo;                  
            jj=i_nodo_vecino(nodos);     
            for nod=1:3
                if i_nodo_vecino(nod)~=ii && i_nodo_vecino(nod)~=jj
                    kk=i_nodo_vecino(nod);
                end
            end               

            % r_i, r_j, r_k: coordenadas de ii,jj,kk
            r_i=[x(ii),y(ii),z(ii)];
            r_j=[x(jj),y(jj),z(jj)]; 
            r_k=[x(kk),y(kk),z(kk)];  

            % es por la contribución doble que está el 1/2? -> por qué hace
            % esto?
            if ii~=jj
                R(i_nodo,jj) = R(i_nodo,jj)+(1/2)*(norm(r_j-r_k,2))/ norm(r_i-r_j,2);
                R(i_nodo,ii) = R(i_nodo,ii) - ((1/2)*(norm(r_j-r_k,2))/norm(r_i-r_j,2));
            end 
        end 

        % area sumada de todos los vecinos (elementos que contienen a un
        % nodo)
        area_total = area_total + area_elemento_i(vecino_de_i);
    end 

    R(i_nodo,:) = (1/area_total)*R(i_nodo,:);  
end  

% LAP = R*F;
LAP = 0;
end


function area=area_triangulo(p1,p2,p3)
    % esta funcion calcula el area de un triangulo definido por su tres
    % vertices
    u=p1-p2;
    v=p1-p3;
    uxv=cross(u,v); %área u-v.
    area=0.5*(norm(uxv,2));
end