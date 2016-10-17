function resultados = driver_searching (Modelo,EGs,fs)
% Assessment of driver location by instant.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUT: 
% Modelo: model of filled atria for 3D representations.
% EGs: signal.
% model: place where rotor occurs {'SR', 'SAF', 'CAF'}
% fs: sampling frequency.
%
%OUTPUT:
% phase: phases for each node.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modelo: estructura de la auricula
% EGs: señal en los nodos [Tiempo x Nodos]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assessment:

% Searching of phase singularities

[singularidades signo_todos angulo_todos]=detecta_singularidades(EGs,Modelo);

resultados.singularidades=singularidades;               %posicion de las singularidades de fase
resultados.signo_todos=signo_todos;                     %sentido de giro de las SP
resultados.angulo_todos=angulo_todos;                   %angulo de giro de las SP

% Seguimiento temporal de los filamentos

[pos_fil sign_fil ang_fil]=filamentos_todos(singularidades,signo_todos,angulo_todos);
[pos_fil_conec_filt sign_fil_conec_filt ang_fil_conec_filt ]=filamentos_vueltas(ang_fil,pos_fil,sign_fil,Modelo,fs);

resultados.pos_rot=pos_fil_conec_filt;                  % posicion de cada rotor a lo largo del tiempo
resultados.sign_rot=sign_fil_conec_filt;                % sentido de giro de cada rotor a lo largo del tiempo
resultados.ang_rot=ang_fil_conec_filt;                  % angulo de cada rotor a lo largo del tiempo
end

function [singularidades signo_todos angulo_todos]=detecta_singularidades(ECGs,Torso);

umbral_ter_cribado=0.15;

Radio1=0.05;
Radio2=0.1;
Radio3=0.15;
Radio4=0.1;
puntos=6;
puntos2=25;

exp_ponderacion=2;

[T n_vertices]=size(ECGs);
singularidades=zeros(T,n_vertices);
signo_todos=zeros(T,n_vertices);
angulo_todos=zeros(T,n_vertices);

angulos=linspace(0,2*pi,puntos+1);
angulos=angulos(1:end-1);
angulos2=linspace(0,2*pi,puntos2+1);
angulos2=angulos2(1:end-1);

vertices=Torso.vertices;
faces=Torso.faces;
dimension=max(vertices(:,1))-min(vertices(:,1));

senyal_radio1=zeros(T,puntos);
senyal_radio2=zeros(T,puntos);
senyal_radio3=zeros(T,puntos);
senyal_radio4=zeros(T,puntos2);

error=zeros(T,n_vertices,3);

for i_p=1:n_vertices
    [a b]=find(faces==i_p);
    vecinos=unique(faces(a,:));
    warning off
    normal = fitNormal(vertices(vecinos,:));
    warning on
    normal=normal./sqrt(sum(normal.^2));
    normal=normal(:)';

    R=sqrt(sum((  ones(n_vertices,1)*vertices(i_p,:) - vertices   ).^2,2));
    puntos_cerca=find(R<= dimension*Radio3*1.5);
    proyeccion=vertices(puntos_cerca,:)- (sum(ones(length(puntos_cerca),1)*normal.*(  vertices(puntos_cerca,:)-ones(length(puntos_cerca),1)*vertices(i_p,:)  ),2))*normal;

    lejos=find(R(puntos_cerca)==max(R(puntos_cerca)));
    vector_x=proyeccion(lejos(1),:)-vertices(i_p,:);
    vector_x=vector_x./sqrt(sum(vector_x.^2));
    vector_y=cross(vector_x,normal);

    puntos_Radio1=ones(puntos,1)*vertices(i_p,:)+dimension.*Radio1.*sin(angulos(:))*vector_x(:)'+dimension.*Radio1.*cos(angulos(:))*vector_y(:)';
    puntos_Radio2=ones(puntos,1)*vertices(i_p,:)+dimension.*Radio2.*sin(angulos(:))*vector_x(:)'+dimension.*Radio2.*cos(angulos(:))*vector_y(:)';
    puntos_Radio3=ones(puntos,1)*vertices(i_p,:)+dimension.*Radio3.*sin(angulos(:))*vector_x(:)'+dimension.*Radio3.*cos(angulos(:))*vector_y(:)';

    proy_inf=[0 0.001 100] - (sum(normal.*(  [0 0.001 100]-vertices(i_p,:)  ),2))*normal;
    vector_x2=proy_inf-vertices(i_p,:);
    vector_x2=vector_x2./sqrt(sum(vector_x2.^2));
    vector_y2=cross(vector_x2,normal);

    puntos_Radio4=ones(puntos2,1)*vertices(i_p,:)+dimension.*Radio4.*sin(angulos2(:))*vector_x2(:)'+dimension.*Radio4.*cos(angulos2(:))*vector_y2(:)';

    for i=1:puntos
       R=sqrt(sum((  ones(length(puntos_cerca),1)*puntos_Radio1(i,:) - proyeccion   ).^2,2));
       pond=1./(R.^exp_ponderacion);
       pond=pond./sum(pond);
       senyal_radio1(:,i)=ECGs(:,puntos_cerca)*pond(:);

       R=sqrt(sum((  ones(length(puntos_cerca),1)*puntos_Radio2(i,:) - proyeccion   ).^2,2));
       pond=1./(R.^exp_ponderacion);
       pond=pond./sum(pond);
       senyal_radio2(:,i)=ECGs(:,puntos_cerca)*pond(:);

       R=sqrt(sum((  ones(length(puntos_cerca),1)*puntos_Radio3(i,:) - proyeccion   ).^2,2));
       pond=1./(R.^exp_ponderacion);
       pond=pond./sum(pond);
       senyal_radio3(:,i)=ECGs(:,puntos_cerca)*pond(:);

%            figure
%            plot(senyal_radio3(:,i))
%            figure
%            plot(ECGs{i_c}(:,puntos_cerca+1).*(ones(T,1)*pond(:)'))
%            pause
%            close all

    end

    for i=1:puntos2
       R=sqrt(sum((  ones(length(puntos_cerca),1)*puntos_Radio4(i,:) - proyeccion   ).^2,2));
       pond=1./(R.^exp_ponderacion);
       pond=pond./sum(pond);
       senyal_radio4(:,i)=ECGs(:,puntos_cerca)*pond(:);
    end

    Fase_radio1=mapadefaseconhilbert100720(senyal_radio1')';
    Fase_radio2=mapadefaseconhilbert100720(senyal_radio2')';
    Fase_radio3=mapadefaseconhilbert100720(senyal_radio3')';
    Fase_radio4=mapadefaseconhilbert100720(senyal_radio4')';

    [vector_sing vector_error vector_signo vector_angulo]=valora_fases(Fase_radio1,Fase_radio2,Fase_radio3,Fase_radio4);         

    singularidades(:,i_p)=vector_sing;
    signo_todos(:,i_p)=vector_signo;
    angulo_todos(:,i_p)=vector_angulo;
    error(:,i_p,:)=vector_error;

%     imagen_torso(Torso,zeros(length(vertices),1));
%     hold on
%     plot3(Torso.vertices(puntos_cerca,1),Torso.vertices(puntos_cerca,2),Torso.vertices(puntos_cerca,3),'.')
%     plot3(Torso.vertices(vecinos,1),Torso.vertices(vecinos,2),Torso.vertices(vecinos,3),'r.','MarkerSize',20)
%     plot3(proyeccion(:,1),proyeccion(:,2),proyeccion(:,3),'k.')
%     plot3(puntos_Radio1(:,1),puntos_Radio1(:,2),puntos_Radio1(:,3),'m.','MarkerSize',20)
%     plot3(puntos_Radio2(:,1),puntos_Radio2(:,2),puntos_Radio2(:,3),'m.','MarkerSize',20)
%     plot3(puntos_Radio3(:,1),puntos_Radio3(:,2),puntos_Radio3(:,3),'m.','MarkerSize',20)
%     quiver3(vertices(i_p,1),vertices(i_p,2),vertices(i_p,3),vector_x(1)*5,vector_x(2)*5,vector_x(3)*5)
%     quiver3(vertices(i_p,1),vertices(i_p,2),vertices(i_p,3),vector_y(1)*5,vector_y(2)*5,vector_y(3)*5)
%     quiver3(vertices(i_p,1),vertices(i_p,2),vertices(i_p,3),normal(1)*5,normal(2)*5,normal(3)*5)
%     axis equal
%     pause
%     close all
end

% eliminamos los que estan cerca en base al error en la linealidad

for i_t=1:T    

    indices=find(singularidades(i_t,:)>0);
    matriz_distancias=ones(length(indices))*1000;
    ind_3=indices;

    for i=1:length(indices)
        for j=i+1:length(indices)
            matriz_distancias(i,j)=sqrt(sum((  Torso.vertices(indices(i),:) - Torso.vertices(indices(j),:)   ).^2)); 
        end
    end

    while ~isempty(find(matriz_distancias(:)<umbral_ter_cribado*dimension))
        [a b]=find(matriz_distancias<umbral_ter_cribado*dimension);
        v_error=[(error(i_t,indices(a(1)),1) < error(i_t,indices(b(1)),1)) ,(error(i_t,indices(a(1)),2) < error(i_t,indices(b(1)),2)),(error(i_t,indices(a(1)),3) < error(i_t,indices(b(1)),3))  ];
        if sum(v_error) > 1  %me quedo con a
            ind_3(b(1))=0;
            matriz_distancias(b(1),:)=1e32;
            matriz_distancias(:,b(1))=1e32;
        else                            %me quedo con b
            ind_3(a(1))=0;
            matriz_distancias(a(1),:)=1e32;
            matriz_distancias(:,a(1))=1e32;
        end
    end 

    a=find(ind_3==0);
    singularidades(i_t,indices(a))=0;
    signo_todos(i_t,indices(a))=0;
    angulo_todos(i_t,indices(a))=0;
end
end

function [pos_fil sign_fil ang_fil]=filamentos_todos(singularidades,signo_todos,angulo_todos);

[T n_vertices]=size(singularidades);
pos_fil=zeros(T,16);
sign_fil=zeros(T,16);
ang_fil=zeros(T,16);

%%%% unimos las sp en filamentos, desde fuera hacia dentro

for i_t=1:T
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% conectividad espacial
    
    p1=find(singularidades(i_t,:)>0);
    signo=signo_todos(i_t,p1);
    angulo=angulo_todos(i_t,p1);

    pos_fil(i_t,1:length(p1))=p1;
    sign_fil(i_t,1:length(p1))=signo;
    ang_fil(i_t,1:length(p1))=angulo;
        
end
end


function [pos_fil_conec_filt sign_fil_conec_filt ang_fil_conec_filt ]=filamentos_vueltas(ang_fil,pos_fil,sign_fil,Torso,fs);

umbral_decision_conec=0.4;
%tiempomin=25;
durmin=0.05;
tiempomin=round(durmin*fs);
vueltas_min=1.2;

dimension=max(Torso.vertices(:,1))-min(Torso.vertices(:,1));

[T nada]=size(pos_fil);

%%%%%% conectividad desde fuera

pos_fil_conec=zeros(T,1000);
sign_fil_conec=zeros(T,1000);
ang_fil_conec=zeros(T,1000);

n_fil=1;
flag_empezar=1;
for i_t=1:T
  
    n_rotores=length(find(pos_fil(i_t,:)>0));
    rotores=find(pos_fil(i_t,:)>0);
    if n_rotores>0
        if flag_empezar
            pos_fil_conec(i_t,n_fil:n_fil+n_rotores-1)=pos_fil(i_t,rotores);
            sign_fil_conec(i_t,n_fil:n_fil+n_rotores-1)=sign_fil(i_t,rotores);
            ang_fil_conec(i_t,n_fil:n_fil+n_rotores-1)=ang_fil(i_t,rotores);
            n_fil=n_fil+n_rotores;
        else
            rotores_ant=find(pos_fil_conec(i_t-1,:)>0);
            matriz_distancias=zeros(length(rotores_ant),n_rotores);
            for i_ra=1:length(rotores_ant)
                for i_r=1:n_rotores
                    matriz_distancias(i_ra,i_r)=sqrt(sum((   Torso.vertices(pos_fil_conec(i_t-1,rotores_ant(i_ra)),:)-Torso.vertices(pos_fil(i_t,rotores(i_r)) ,:)).^2));
                end
            end
            
            n_asig=min(n_rotores,length(rotores_ant));
            vector_faltan=zeros(n_rotores,1);
            
            for i_a=1:n_asig
                [rot_ant,rot_ac]=find(matriz_distancias==min(min(matriz_distancias)));
                rot_ant=rot_ant(1);
                rot_ac=rot_ac(1);
                if matriz_distancias(rot_ant,rot_ac)<umbral_decision_conec*dimension
                    pos_fil_conec(i_t,rotores_ant(rot_ant))=pos_fil(i_t,rotores(rot_ac));
                    sign_fil_conec(i_t,rotores_ant(rot_ant))=sign_fil(i_t,rotores(rot_ac));
                    ang_fil_conec(i_t,rotores_ant(rot_ant))=ang_fil(i_t,rotores(rot_ac));
                    matriz_distancias(rot_ant,:)=1e12*ones(1,n_rotores);
                    matriz_distancias(:,rot_ac)=1e12*ones(length(rotores_ant),1);
                    vector_faltan(rot_ac)=1;
                else
                    break
                end
            end
            faltan=find(vector_faltan==0);
            if ~isempty(faltan)
                pos_fil_conec(i_t,n_fil:n_fil+length(faltan)-1,:)=pos_fil(i_t,rotores(faltan));
                sign_fil_conec(i_t,n_fil:n_fil+length(faltan)-1,:)=sign_fil(i_t,rotores(faltan));
                ang_fil_conec(i_t,n_fil:n_fil+length(faltan)-1,:)=ang_fil(i_t,rotores(faltan));
                n_fil=n_fil+length(faltan);
            end
        end
        flag_empezar=0;
    else
        flag_empezar=1;
    end
end

%%%%%% mirar si da una vuelta completa

n_filamentos=sum( sum(pos_fil_conec,1)>0);
fil_buenos=zeros(n_filamentos,1);

for i=1:n_filamentos
    
    ts_fuera=find(pos_fil_conec(:,i)>0);
    vector_angulos=unwrap(ang_fil_conec(ts_fuera,i));
%          
%     figure
%     plot(vector_angulos)
%     pause
%     close all
    
    if length(vector_angulos)>tiempomin
        if max(vector_angulos)-min(vector_angulos) > 2*pi*vueltas_min
            fil_buenos(i)=1;
        end
   end
end

fil_buenos=find(fil_buenos==1);
pos_fil_conec_filt = pos_fil_conec(:,fil_buenos);
ang_fil_conec_filt = ang_fil_conec(:,fil_buenos);
sign_fil_conec_filt = sign_fil_conec(:,fil_buenos);
end


function n = fitNormal(data, show_graph)
%FITNORMAL - Fit a plane to the set of coordinates
%
%For a passed list of points in (x,y,z) cartesian coordinates,
%find the plane that best fits the data, the unit vector
%normal to that plane with an initial point at the average
%of the x, y, and z values.
%
% :param data: Matrix composed of of N sets of (x,y,z) coordinates
%              with dimensions Nx3
% :type data: Nx3 matrix
%
% :param show_graph: Option to display plot the result (default false)
% :type show_graph: logical
%
% :return n: Unit vector that is normal to the fit plane
% :type n: 3x1 vector
	
	if nargin == 1
		show_graph = false;
	end
	
	for i = 1:3
		X = data;
		X(:,i) = 1;
		
		X_m = X' * X;
		if det(X_m) == 0
			can_solve(i) = 0;
			continue
		end
		can_solve(i) = 1;
		
		% Construct and normalize the normal vector
		coeff = (X_m)^-1 * X' * data(:,i);
		c_neg = -coeff;
		c_neg(i) = 1;
		coeff(i) = 1;
		n(:,i) = c_neg / norm(coeff);
		
	end
	
	if sum(can_solve) == 0
		error('Planar fit to the data caused a singular matrix.')
		return
	end
	
	% Calculating residuals for each fit
	center = mean(data);
	off_center = [data(:,1)-center(1) data(:,2)-center(2) data(:,3)-center(3)];
	for i = 1:3
		if can_solve(i) == 0
			residual_sum(i) = NaN;
			continue
		end
		
		residuals = off_center * n(:,i);
		residual_sum(i) = sum(residuals .* residuals);
		
	end
	
	% Find the lowest residual index
	best_fit = find(residual_sum == min(residual_sum));
	
	% Possible that equal mins so just use the first index found
	n = n(:,best_fit(1));
	
	if ~show_graph
		return
	end
	
	range = max(max(data) - min(data)) / 2;
	mid_pt = (max(data) - min(data)) / 2 + min(data);
	xlim = [-1 1]*range + mid_pt(1);
	ylim = [-1 1]*range + mid_pt(2);
	zlim = [-1 1]*range + mid_pt(3);

	L=plot3(data(:,1),data(:,2),data(:,3),'ro','Markerfacecolor','r'); % Plot the original data points
	hold on;
	set(get(L, 'Parent'),'DataAspectRatio',[1 1 1],'XLim',xlim,'YLim',ylim,'ZLim',zlim);
	
	norm_data = [mean(data); mean(data) + (n' * range)];
	
	% Plot the original data points
	L=plot3(norm_data(:,1),norm_data(:,2),norm_data(:,3),'b-','LineWidth',3);
	set(get(get(L,'parent'),'XLabel'),'String','x','FontSize',14,'FontWeight','bold')
	set(get(get(L,'parent'),'YLabel'),'String','y','FontSize',14,'FontWeight','bold')
	set(get(get(L,'parent'),'ZLabel'),'String','z','FontSize',14,'FontWeight','bold')
	title(sprintf('Normal Vector: <%0.3f, %0.3f, %0.3f>',n),'FontWeight','bold','FontSize',14)
	grid on;
	axis square;
	hold off;
end

function [F M]=mapadefaseconhilbert100720(A)

%warning off
  if length(size(A))==3
              [nt,nx,ny]=size(A);
              F=zeros(size(A));
              M=zeros(size(A));
             for i=1:nx
                 for j=1:ny
                    s=squeeze(A(:,i,j));
                    s=s-mean(s);
                    h=hilbert(s);
                    ha=imag(h);
                    F(:,i,j)=-atan2(ha,s);
                    M(:,i,j)=abs(h);
                 end
             end

            
  end
  
  
  if length(size(A))==2
        [np,nt]=size(A);
        F=zeros(size(A));
        M=zeros(size(A));
        
        for i=1:np
            s=A(i,:)-mean(A(i,:));
            h=hilbert(s);
            ha=imag(h);
            F(i,:)=-atan2(ha,s); 
            M(i,:)=abs(h);
        end
        
  end
end


function [vector_sing vector_error vector_signo vector_angulo]=valora_fases(Fase_radio1,Fase_radio2,Fase_radio3,Fase_radio4);

umbral_primer_cribado=0.4;
umbral_seg_cribado=0.4;

vector_sing=zeros(size(Fase_radio1,1),1);
vector_signo=zeros(size(Fase_radio1,1),1);
vector_angulo=zeros(size(Fase_radio1,1),1);
vector_error=zeros(size(Fase_radio1,1),3);

difFase1=abs(diff([Fase_radio1 Fase_radio1(:,1)],1,2));
difFase2=abs(diff([Fase_radio2 Fase_radio2(:,1)],1,2));
difFase3=abs(diff([Fase_radio3 Fase_radio3(:,1)],1,2));

ind_1=find(sum(difFase1>=umbral_primer_cribado*2*pi,2));

if ~isempty(ind_1)
    
    v1=sum( difFase1(ind_1,:)>=umbral_seg_cribado*2*pi,2);
    v2=sum(difFase2(ind_1,:)>=umbral_seg_cribado*2*pi,2);
    v3=sum(difFase3(ind_1,:)>=umbral_seg_cribado*2*pi,2);
    vdec=[v1 v2 v3];
    
    ind_2=find(sum(vdec==1,2)>1);
    
%     for i=1:length(ind_1)
%         figure
%         hold on
%         plot(Fase_radio1(ind_1(i),:))
%         plot(Fase_radio2(ind_1(i),:),'r')
%         plot(Fase_radio3(ind_1(i),:),'g')
%         pause
%         close all
%     end
        
   
    if ~isempty(ind_2)
        
%         vector_sing(ind_1)=1;
        vector_sing(ind_1(ind_2))=1;
        
        for i=1:length(ind_2)
            fase=Fase_radio1(ind_1(ind_2(i)),:);
            MF_dif=fase-[fase(end) fase(1:end-1)];
            ind=find(abs(MF_dif)==max(abs(MF_dif)));
            MF_lineal=[fase(ind:end) fase(1:ind-1)];
            X=[ones(length(fase),1) (1:length(fase))'];
            Y=MF_lineal';
            M=(X'*X)\(X'*Y);
            vector_error(ind_1(ind_2(i)),1)=sum((Y-X*M).^2);
            signo1=sign(M(2));
            
            fase=Fase_radio2(ind_1(ind_2(i)),:);
            MF_dif=fase-[fase(end) fase(1:end-1)];
            ind=find(abs(MF_dif)==max(abs(MF_dif)));
            MF_lineal=[fase(ind:end) fase(1:ind-1)];
            X=[ones(length(fase),1) (1:length(fase))'];
            Y=MF_lineal';
            M=(X'*X)\(X'*Y);
            vector_error(ind_1(ind_2(i)),2)=sum((Y-X*M).^2);
            signo2=sign(M(2));
            
            fase=Fase_radio3(ind_1(ind_2(i)),:);
            MF_dif=fase-[fase(end) fase(1:end-1)];
            ind=find(abs(MF_dif)==max(abs(MF_dif)));
            MF_lineal=[fase(ind:end) fase(1:ind-1)];
            X=[ones(length(fase),1) (1:length(fase))'];
            Y=MF_lineal';
            M=(X'*X)\(X'*Y);
            vector_error(ind_1(ind_2(i)),3)=sum((Y-X*M).^2);
            signo3=sign(M(2));
            
            vector_signo(ind_1(ind_2(i)))=sign(signo1+signo2+signo3);
            
            fase=Fase_radio4(ind_1(ind_2(i)),:);
            MF_dif=fase-[fase(end) fase(1:end-1)];
            indice_angulo=find(abs(MF_dif)==max(abs(MF_dif)));            
            vector_angulo(ind_1(ind_2(i)))=indice_angulo/length(fase)*2*pi;
        end

%         for i=1:length(ind_2)
%             figure
%             hold on
%             plot(Fase_radio1(ind_1(ind_2(i)),:))
%             plot(Fase_radio2(ind_1(ind_2(i)),:),'r')
%             plot(Fase_radio3(ind_1(ind_2(i)),:),'g')
%             pause
%             close all
%         end
        
    end
end
end



%%%%%%%%%%%%%%
%%%%%%%%%%%%%%
%%% Miguel Rodrigo Bort, ITACA, UPV, 2015
%%%%%%%%%%%%%%
%%%%%%%%%%%%%%