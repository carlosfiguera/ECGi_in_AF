function [vecinos] = compile (elementos)
% esta funcion se encarga de...
% dada una matriz cuyas filas los indices de los nodos que constituyen
% cada elemento a cada nodo le asocia el numero de elemento que lo 
% continen y el indice de dichos elemento
% se entiendo que la numeracion de los nodos va desde 1:1:maximo_indice
[num_elementos,aux]=size(elementos);
num_nodos=max(max(elementos));

%%% nota supondremos que tenemos un maximo de 14 amigos de un nodo 
%%% si fuese necesario se puede cambiar facilmente
vecinos=zeros(num_nodos,15);

for elemento=1:num_elementos
n_e=elementos(elemento,:); 
    for kk=1:3
        vecinos(n_e(kk),1)=vecinos(n_e(kk),1)+1;    %numero de faces en las que está contenido el nodo.
        vecinos(n_e(kk),vecinos(n_e(kk),1)+1)=elemento; %faces que contienen al nodo.
    end
end