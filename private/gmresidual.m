function x_hat = gmresidual (A, y)
% GMRES reconstruction method.
%
% This routine, by victor.suarez.gutierrez@urjc.es


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUT:
% A: transfer matrix.
% y: torso potentials.

%OUTPUT:
% x_hat: estimated epicardial potentials.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


A2 = A'*A;
tol = 0.0001;
iter = 30;
for j=1:length(y(1,:))
    fprintf('iteracion nº %i\n', j);
    b2 = A'*y(:,j);
    [x,~,~,~] = gmres (A2,b2,[],tol,iter);
    x_hat(:,j)=x;
end








