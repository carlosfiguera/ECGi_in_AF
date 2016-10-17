function [AA, L, LL] = precompute_matrices (A, order, atrial_model)
% Function which applies narrow-band filter to x and calculates variables for 
% minimizing computing time in calculation of inverse matrix A.
%
% This routine by Víctor Suárez Gutiérrez (victor.suarez.gutierrez@urjc.es)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUT:
% A: transfer matrix.
% order: Tikhonov order (0, 1 or 2).
% atrial_model: atrial model.
%
%OUTPUT:
% AA: A'*A
% L: matrix which takes part in regularization term of Tikhonov approach.
% LL: L'*L
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% minimizing computing time:
    if ~isnan(order)
        [AA, L, LL] = variables_tik (A, order, atrial_model);
    else
        AA = [];
        L = [];
        LL = [];
    end
end


function [AA, L, LL] = variables_tik (A, order, atrial_model)
%Routine for calculating variables used to Tikhonov approaches: A'*A y L'*L.
%
% This routine by Víctor Suárez Gutiérrez (victor.suarez.gutierrez@urjc.es)

%%
switch order
    case 0,     
        [~,T] = size(A);
        L = eye(T); % identity matrix.
    case 1,
        xv  = atrial_model.vertices(:,1);
        yv  = atrial_model.vertices(:,2);
        zv  = atrial_model.vertices(:,3);
        tri = atrial_model.faces; 
        L = gradiente_mat(xv,yv,zv,tri,A);
    case 2,
        xv  = atrial_model.vertices(:,1);
        yv  = atrial_model.vertices(:,2);
        zv  = atrial_model.vertices(:,3);
        tri = atrial_model.faces; 
        L = laplaciano_mat(xv,yv,zv,tri,A);
end

AA = A'*A;
LL = L'*L;

end


function [R,GRAD] = gradiente_mat(x,y,z,elements,F)
% Assessment of gradient.

% Inicialization:
[~,N] = size(F);
R     = zeros(N,N);

[num_ele, ~]   = size(elements);
num_nodes       = size(x);
area_element_i = zeros(num_ele,1);
neighbour = compile (elements);   

% Compute area of each triangle.
for element=1:num_ele
    n_e = elements(element,:);    % nodes places over vertex.
    P1=[x(n_e(1)),y(n_e(1)),z(n_e(1))];     % coordinates of each node.
    P2=[x(n_e(2)),y(n_e(2)),z(n_e(2))];
    P3=[x(n_e(3)),y(n_e(3)),z(n_e(3))];
    area_element_i(element) = triangle_area (P1,P2,P3); % triangle area of these thre points.
end
    
% Create matrix B
% by node
for i_node=1:num_nodes
    area_total=0; 
    num_neighbour = neighbour (i_node,1);

    % each neighbour of a node
    for i_neighbour=1:num_neighbour;

        neighbour_i = neighbour (i_node,i_neighbour+1);
        i_neighbour_node = elements (neighbour_i,:); % nodes do with each face 'neighbour_i'.

        for nodes=1:3 
            % ii, jj, kk interesting: index for different nodes of element ii
            ii=i_node;                  
            jj=i_neighbour_node(nodes);     
            for nod=1:3
                if i_neighbour_node(nod)~=ii && i_neighbour_node(nod)~=jj
                    kk=i_neighbour_node(nod);
                end
            end               

            % r_i, r_j, r_k: coordinates of ii,jj,kk
            r_i=[x(ii),y(ii),z(ii)];
            r_j=[x(jj),y(jj),z(jj)]; 
            r_k=[x(kk),y(kk),z(kk)];  

            if ii~=jj
                R(i_node,jj) = R(i_node,jj)+(1/2)*(norm(r_j-r_k,2))/ norm(r_i-r_j,2);
                R(i_node,ii) = R(i_node,ii) - ((1/2)*(norm(r_j-r_k,2))/norm(r_i-r_j,2));
            end 
        end 

        area_total = area_total + area_element_i(neighbour_i);
    end 
end  
GRAD = 0;
end


function [R,LAP]=laplaciano_mat(x,y,z,elements,F)
% Assessment of laplacian.

[~,N] = size(F);
R = zeros(N,N);

[num_ele aux] = size(elements);
num_nodes = size(x);
area_element_i = zeros(num_ele,1);
neighbour = compile (elements);

% Compute area of each triangle.
for element=1:num_ele
    n_e=elements(element,:);
    P1=[x(n_e(1)),y(n_e(1)),z(n_e(1))];
    P2=[x(n_e(2)),y(n_e(2)),z(n_e(2))];
    P3=[x(n_e(3)),y(n_e(3)),z(n_e(3))];
    area_element_i(element)=triangle_area(P1,P2,P3);
end
    
% Create matrix B
% by node
for i_node=1:num_nodes
    area_total=0; 
    num_neighbour = neighbour (i_node,1);

    % each neighbour of a node
    for i_neighbour=1:num_neighbour;

        neighbour_i = neighbour (i_node,i_neighbour+1);
        i_neighbour_node = elements (neighbour_i,:);

        for nodes=1:3 

            % ii, jj, kk interesting: index for different nodes of element ii
            ii=i_node;                  
            jj= i_neighbour_node (nodes);     
            for nod=1:3
                if i_neighbour_node(nod)~=ii && i_neighbour_node(nod)~=jj
                    kk=i_neighbour_node(nod);
                end
            end               

            % r_i, r_j, r_k: coordinates of ii,jj,kk
            r_i=[x(ii),y(ii),z(ii)];
            r_j=[x(jj),y(jj),z(jj)]; 
            r_k=[x(kk),y(kk),z(kk)];  

            if ii~=jj
                R(i_node,jj) = R(i_node,jj)+(1/2)*(norm(r_j-r_k,2))/ norm(r_i-r_j,2);
                R(i_node,ii) = R(i_node,ii) - ((1/2)*(norm(r_j-r_k,2))/norm(r_i-r_j,2));
            end 
        end 

        area_total = area_total + area_element_i(neighbour_i);
    end 

    R(i_node,:) = (1/area_total)*R(i_node,:);  
end  

LAP = 0;
end


function area = triangle_area (p1,p2,p3)
% This function calculates the area of each face (triangle) from three points.
% p1, p2 and p3.

    u=p1-p2;
    v=p1-p3;
    uxv=cross(u,v); % area u-v.
    area=0.5*(norm(uxv,2));
end