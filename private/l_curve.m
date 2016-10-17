function [reg_corner, rho, eta, reg_param] = l_curve (U, sm, sk, b, method, SNR)
%L_CURVE Plot the L-curve and find its "corner".
%
% This routine by Víctor Suárez Gutiérrez (victor.suarez.gutierrez@urjc.es)
% based on Christian Hansen´s algorithm (DTU Compute, October 27, 2010).
%
%
%%
% [reg_corner,rho,eta,reg_param] =
%                  l_curve(U,s,b,method)
%                  l_curve(U,sm,b,method)  ,  sm = [sigma,mu]
%                  l_curve(U,s,b,method,L,V)
%
% Plots the L-shaped curve of eta, the solution norm || x || or
% semi-norm || L x ||, as a function of rho, the residual norm
% || A x - b ||, for the following methods:
%    method = 'Tikh'  : Tikhonov regularization   (solid line )
%    method = 'grst'  : Greensite regularization   (solid line )
%    method = 'tsvd'  : truncated SVD or GSVD     (o markers  )
%    method = 'dsvd'  : damped SVD or GSVD        (dotted line)
% The corresponding reg. parameters are returned in reg_param.  If no
% method is specified then 'Tikh' is default.  For other methods use plot_lc.
%
% Note that 'Tikh', 'tsvd' and 'dsvd' require either U and s (standard-
% form regularization) computed by the function csvd, or U and sm (general-
% form regularization) computed by the function cgsvd, while 'mtvsd'
% requires U and s as well as L and V computed by the function csvd.
%
% If any output arguments are specified, then the corner of the L-curve
% is identified and the corresponding reg. parameter reg_corner is
% returned.  Use routine l_corner if an upper bound on eta is required.


if (nargin==3)
    method='Tikh'; 
end  % Tikhonov reg. is default.
npoints = 500;  % Number of points on the L-curve for Tikh and dsvd.
smin_ratio = 16*eps;  % Smallest regularization parameter.

% Initialization.
[m,n] = size(U); [p,ps] = size(sm);
if (nargout > 0)
    locate = 1; 
else
    locate = 0; 
end
beta = U'*b; 
beta2 = norm(b)^2 - norm(beta)^2;


if (strncmp(method,'Tikh',4) || strncmp(method,'tikh',4) || strncmp(method,'grst',4))

    if (ps==1)
      s = sm; beta = beta(1:p);
      s1 = sk; 
    else
      s = sm(p:-1:1,1)./sm(p:-1:1,2); beta = beta(p:-1:1);
      s1 = sk(p:-1:1,1)./sk(p:-1:1,2);
    end
    xi = beta(1:p)./s1;
    xi( isinf(xi) ) = 0;
    
    
  eta = zeros(npoints,1); 
  rho = eta; reg_param = eta; 
  s2 = s.^2;
  if SNR~=40
    reg_param(npoints) = 1e-15;   
  else
    reg_param(npoints) = 1e-20;
  end
  ratio = (s1(1)/reg_param(npoints))^(1/(npoints-1));
  for i=npoints-1:-1:1, 
      reg_param(i) = ratio*reg_param(i+1); 
  end
  for i=1:npoints
    f = s2./(s2 + reg_param(i)^2);
    eta(i) = norm(f.*xi);
    rho(i) = norm((1-f).*beta(1:p));
  end
  if (m > n && beta2 > 0)
      rho = sqrt(rho.^2 + beta2); 
  end
  marker = '-'; txt = 'Tikh.';

elseif (strncmp(method,'tsvd',4) || strncmp(method,'tgsv',4))
    if (ps==1)
      s = sm; beta = beta(1:p);
      s1 = sk; 
    else
      s = sm(p:-1:1,1)./sm(p:-1:1,2); beta = beta(p:-1:1);
      s1 = sk(p:-1:1,1)./sk(p:-1:1,2);
    end
    xi = beta(1:p)./s1;
    xi( isinf(xi) ) = 0;
    
    
  eta = zeros(p,1); rho = eta;
  eta(1) = abs(xi(1))^2;
  for k=2:p
      eta(k) = eta(k-1) + abs(xi(k))^2;
  end
  eta = sqrt(eta);
  if (m > n)
    if (beta2 > 0)
        rho(p) = beta2; 
    else
        rho(p) = eps^2;
    end
  else
    rho(p) = eps^2;
  end
  for k=p-1:-1:1
      rho(k) = rho(k+1) + abs(beta(k+1))^2; 
  end
  rho = sqrt(rho);
  reg_param = (1:p)'; marker = 'o';
  if (ps==1)
    U = U(:,1:p); txt = 'TSVD';
  else
    U = U(:,1:p); txt = 'TGSVD';
  end

elseif (strncmp(method,'dsvd',4) || strncmp(method,'dgsv',4))
    if (ps==1)
      s = sm; 
      s1 = sk; 
    else
      [s,s_pos] = sort(sm(:,1),'descend');
      s1 = sm(s_pos,2);
    end
    xi = beta(1:p)./s;
    xi( isinf(xi) ) = 0;
    
    
  eta = zeros(npoints,1); 
  rho = eta; reg_param = eta;
  reg_param(npoints) = 1e-15;   
  ratio = (s(1)/reg_param(npoints))^(1/(npoints-1));
  for i=npoints-1:-1:1
      reg_param(i) = ratio*reg_param(i+1); 
  end
  for i=1:npoints
    f = s./(s + reg_param(i).*s1);
    eta(i) = norm(f.*xi);
    rho(i) = norm((1-f).*beta(1:p));
  end
  if (m > n & beta2 > 0)
      rho = sqrt(rho.^2 + beta2); 
  end
  marker = ':';
  if (ps==1)
      txt = 'DSVD'; 
  else
      txt = 'DGSVD'; 
  end

else
  error('Illegal method')
end

% Locate the "corner" of the L-curve, if required.
if (locate)
    [reg_corner, rho_c, eta_c, selection] = l_corner (rho, eta, reg_param, U, sm, sk, b, method);
end

if 0    % change 0 for 1 for L-curve representation.
    % Make plot.
    plot_lc(rho,eta,marker,ps,reg_param);
    if locate
      ax = axis;
      HoldState = ishold; hold on;
      loglog([min(rho)/100,rho_c],[eta_c,eta_c],':r',...
             [rho_c,rho_c],[min(eta)/100,eta_c],':r')
      title(['L-curve, ',txt,' corner at ',num2str(reg_corner)]);
      axis(ax)
      if (~HoldState), hold off; end
    end
end
end


function [reg_c, rho_c, eta_c, selection] = l_corner (rho, eta, reg_param, U, s, s1, b, method)
%L_CORNER Locate the "corner" of the L-curve.

% This routine by Víctor Suárez Gutiérrez (victor.suarez.gutierrez@urjc.es)
% based on Christian Hansen´s algorithm (DTU Compute, October 27, 2010).


%%
% [reg_c,rho_c,eta_c] =
%        l_corner(rho,eta,reg_param)
%        l_corner(rho,eta,reg_param,U,s,b,method,M)
%        l_corner(rho,eta,reg_param,U,sm,b,method,M) ,  sm = [sigma,mu]
%
% Locates the "corner" of the L-curve in log-log scale.
%
% It is assumed that corresponding values of || A x - b ||, || L x ||,
% and the regularization parameter are stored in the arrays rho, eta,
% and reg_param, respectively (such as the output from routine l_curve).
%
% If nargin = 3, then no particular method is assumed, and if
% nargin = 2 then it is issumed that reg_param = 1:length(rho).
%
% If nargin >= 6, then the following methods are allowed:
%    method = 'Tikh'  : Tikhonov regularization
%    method = 'grst'  : Greensite regularization   (solid line )
%    method = 'tsvd'  : truncated SVD or GSVD
%    method = 'dsvd'  : damped SVD or GSVD
%    method = 'mtsvd' : modified TSVD,
% and if no method is specified, 'Tikh' is default.  If the Spline Toolbox
% is not available, then only 'Tikh' and 'dsvd' can be used.
%
% An eighth argument M specifies an upper bound for eta, below which
% the corner should be found.


% Ensure that rho and eta are column vectors.
rho = rho(:); 
eta = eta(:);

% Set default regularization method.
if (nargin <= 3)
  method = 'none';
  if (nargin==2)
      reg_param = (1:length(rho))'; 
  end
else
  if (nargin==6)
      method = 'Tikh'; 
  end
end

% Set this logical variable to 1 (true) if the corner algorithm
% should always be used, even if the Spline Toolbox is available.
alwayscorner = 0;

% Set threshold for skipping very small singular values in the
% analysis of a discrete L-curve.
s_thr = eps;  % Neglect singular values less than s_thr.

% Set default parameters for treatment of discrete L-curve.
deg   = 2;  % Degree of local smooting polynomial.
q     = 2;  % Half-width of local smoothing interval.
order = 4;  % Order of fitting 2-D spline curve.

% Initialization.

if (nargin > 3)
  [p,ps] = size(s); [m,n] = size(U);
  beta = U'*b; b0 = b - U*beta;
  if (ps==2)
    s = s(p:-1:1,1)./s(p:-1:1,2);
    s1 = s1(p:-1:1,1)./s1(p:-1:1,2);
    beta = beta(p:-1:1);
  end
  xi = beta./s1;
  if (m>n)  % Take of the least-squares residual.
      beta = [beta;norm(b0)];
  end
end

if (strncmp(method,'Tikh',4) || strncmp(method,'tikh',4)) || strncmp(method,'grst',4)
  % The L-curve is differentiable; computation of curvature in
  % log-log scale is easy.
    x = log10(rho);
    y = log10(eta);
    dx = gradient(x,2); %uses h (x,h) as the spacing between points in each direction.
    dy = gradient(y,2); %uses h (x,h) as the spacing between points in each direction.

    %optimal alpha vector:
    slope = (dy./dx);
    [MM ll] = max(slope);
    ddx = gradient(dx,2);
    ddy = gradient(dy,2);
    curva =(dx.*ddy-ddx.*dy)./((dx.^2+dy.^2).^(3/2));
    maxcurva  = find(abs(curva) == max(abs(curva)), 1, 'first');
    if maxcurva<=0.5*length(rho)
        curva2=curva(maxcurva+0.1*length(rho):end);
        curvainv=curva2(end:-1:1);
        maxcurva  = find(abs(curva) == max(abs(curva)), 1, 'first');
        maxcurva2  = length(rho)-find(abs(curvainv) == max(abs(curvainv)), 1, 'first');
    else
        countt=0;
        while maxcurva>0.3*length(rho) && countt<=20
            curva2=curva(1:maxcurva-round(countt*0.05*length(rho)));
            curvainv=curva2(end:-1:1);
            maxcurva2  = maxcurva-round(countt*0.05*length(rho));
            maxcurva  = maxcurva2-find(abs(curvainv) == max(abs(curvainv)), 1, 'first');
            countt = countt+1;
        end
    end
    grados=atan(dy./dx)*(180/pi);
    grados=grados(maxcurva:maxcurva2);
    diferencia=diff(grados);
    [mm, ~] = max(abs(diferencia(1:round(length(diferencia)/2))));
    optimal = find(abs(diferencia) <= mm*0.33, 1, 'first');
    if isempty(optimal)
        [~, optimal] = min(diferencia);
    end
    reg_c = reg_param(optimal+maxcurva);
    f = (s.^2)./(s.^2 + reg_c^2);
    eta_c = norm(f.*xi);
    rho_c = norm((1-f).*beta(1:length(f)));
    selection=optimal+maxcurva;
        
elseif (strncmp(method,'tsvd',4) || strncmp(method,'tgsv',4) || ...
        strncmp(method,'mtsv',4) || strncmp(method,'none',4))

  % Use the adaptive pruning algorithm to find the corner, if the
  % Spline Toolbox is not available.
  if ~exist('splines','dir') || alwayscorner
    %error('The Spline Toolbox in not available so l_corner cannot be used')
    reg_c = corner(rho,eta);
    rho_c = rho(reg_c);
    eta_c = eta(reg_c);
    return
  end

  % Otherwise use local smoothing followed by fitting a 2-D spline curve
  % to the smoothed discrete L-curve. Restrict the analysis of the L-curve
  % according to s_thr.
  if (nargin > 3)    
    index = find(s > s_thr);
    rho = rho(index); 
    eta = eta(index); 
    reg_param = reg_param(index);
  end

  % Convert to logarithms.
  lr = length(rho);
  lrho = log(rho); 
  leta = log(eta); 
  slrho = lrho; 
  sleta = leta;

  % For all interior points k = q+1:length(rho)-q-1 on the discrete
  % L-curve, perform local smoothing with a polynomial of degree deg
  % to the points k-q:k+q.

  v = (-q:q)'; 
  A = zeros(2*q+1,deg+1); 
  A(:,1) = ones(length(v),1);
  for j = 2:deg+1
      A(:,j) = A(:,j-1).*v; 
  end
  for k = q+1:lr-q-1
    cr = A\lrho(k+v); slrho(k) = cr(1);
    ce = A\leta(k+v); sleta(k) = ce(1);
  end

  % Fit a 2-D spline curve to the smoothed discrete L-curve.
  sp = spmak((1:lr+order),[slrho';sleta']);
  pp = ppbrk(sp2pp(sp),[4,lr+1]);

  % Extract abscissa and ordinate splines and differentiate them.
  % Compute as many function values as default in spleval.
  P     = spleval(pp);  dpp   = fnder(pp);
  D     = spleval(dpp); ddpp  = fnder(pp,2);
  DD    = spleval(ddpp);
  ppx   = P(1,:);       ppy   = P(2,:);
  dppx  = D(1,:);       dppy  = D(2,:);
  ddppx = DD(1,:);      ddppy = DD(2,:);

  % Compute the corner of the discretized .spline curve via max. curvature.
  % No need to refine this corner, since the final regularization
  % parameter is discrete anyway.
  % Define curvature = 0 where both dppx and dppy are zero.

  k1    = dppx.*ddppy - ddppx.*dppy;
  k2    = (dppx.^2 + dppy.^2).^(1.5);
  I_nz  = find(k2 ~= 0);
  kappa = zeros(1,length(dppx));
  kappa(I_nz) = -k1(I_nz)./k2(I_nz);
  [kmax , ~] = max(kappa); 
  grados = atan(diff(ppy)./diff(ppx))*(180/pi);
  grados(grados>-85) = 0;
  kappa1 = kappa;
  kappa1(kappa1~=0) = 1;
  kappa1(1:10)=0;
  resultado = kappa1(1:end-1).*grados;
  res = find(resultado~=0);
  x_corner = ppx(res(1)); y_corner = ppy(res(1));
    
  % Locate the point on the discrete L-curve which is closest to the
  % corner of the spline curve.  Prefer a point below and to the
  % left of the corner.  If the curvature is negative everywhere,
  % then define the leftmost point of the L-curve as the corner.
  if (kmax < 0)
    reg_c = reg_param(lr); 
    rho_c = rho(lr); 
    eta_c = eta(lr);
  else
    index = find(lrho < x_corner & leta < y_corner);
    if ~isempty(index)
      [~,rpi] = min((lrho(index)-x_corner).^2 + (leta(index)-y_corner).^2);
      rpi = index(rpi);
    else
      [~,rpi] = min((lrho-x_corner).^2 + (leta-y_corner).^2);
    end
    reg_c = reg_param(rpi); 
    rho_c = rho(rpi); 
    eta_c = eta(rpi);
  end
  selection=0;

elseif (strncmp(method,'dsvd',4) || strncmp(method,'dgsv',4))

  % Use the adaptive pruning algorithm to find the corner, if the
  % Spline Toolbox is not available.ç
  if ~exist('splines','dir') || alwayscorner
    %error('The Spline Toolbox in not available so l_corner cannot be used')
    reg_c = corner(rho,eta);
    rho_c = rho(reg_c);
    eta_c = eta(reg_c);
    return
  end

  % Otherwise use local smoothing followed by fitting a 2-D spline curve
  % to the smoothed discrete L-curve. Restrict the analysis of the L-curve
  % according to s_thr.
%   if (nargin > 3)    
%     index = find(s > s_thr);
%     rho = rho(index); 
%     eta = eta(index); 
%     reg_param = reg_param(index);
%   end

  % Convert to logarithms.
  lr = length(rho);
  lrho = log(rho); 
  leta = log(eta); 
  slrho = lrho; 
  sleta = leta;

  % For all interior points k = q+1:length(rho)-q-1 on the discrete
  % L-curve, perform local smoothing with a polynomial of degree deg
  % to the points k-q:k+q.

  v = (-q:q)'; 
  A = zeros(2*q+1,deg+1); 
  A(:,1) = ones(length(v),1);
  for j = 2:deg+1
      A(:,j) = A(:,j-1).*v; 
  end
  for k = q+1:lr-q-1
    cr = A\lrho(k+v); slrho(k) = cr(1);
    ce = A\leta(k+v); sleta(k) = ce(1);
  end

  % Fit a 2-D spline curve to the smoothed discrete L-curve.
  sp = spmak((1:lr+order),[slrho';sleta']);
  pp = ppbrk(sp2pp(sp),[4,lr+1]);

  % Extract abscissa and ordinate splines and differentiate them.
  % Compute as many function values as default in spleval.
  P     = spleval(pp);  dpp   = fnder(pp);
  D     = spleval(dpp); ddpp  = fnder(pp,2);
  DD    = spleval(ddpp);
  ppx   = P(1,:);       ppy   = P(2,:);
  dppx  = D(1,:);       dppy  = D(2,:);
  ddppx = DD(1,:);      ddppy = DD(2,:);

  % Compute the corner of the discretized .spline curve via max. curvature.
  % No need to refine this corner, since the final regularization
  % parameter is discrete anyway.
  % Define curvature = 0 where both dppx and dppy are zero.
% keyboard;
  k1    = dppx.*ddppy - ddppx.*dppy;
  k2    = (dppx.^2 + dppy.^2).^(1.5);
  I_nz  = find(k2 ~= 0);
  kappa = zeros(1,length(dppx));
  kappa(I_nz) = -k1(I_nz)./k2(I_nz);
  [kmax , ~] = max(kappa); 
  grados = atan(diff(ppy)./diff(ppx))*(180/pi);
  grados(grados>-85) = 0;
  if sum(grados)==0 % in order 0, avoiding empty vector error.
    grados = atan(diff(ppy)./diff(ppx))*(180/pi);
    grados(grados>0) = 0;
    [~, bbb] = min(grados);
    vector = zeros(1,length(grados));
    vector(bbb) = grados(bbb); 
  else
    grados(grados>-85) = 0;
  end
  kappa1 = kappa;
  kappa1(kappa1~=0) = 1;
  kappa1(1:10)=0;
  resultado = kappa1(1:end-1).*grados;
  res = find(resultado~=0);
  x_corner = ppx(res(1)); y_corner = ppy(res(1));
    
  % Locate the point on the discrete L-curve which is closest to the
  % corner of the spline curve.  Prefer a point below and to the
  % left of the corner.  If the curvature is negative everywhere,
  % then define the leftmost point of the L-curve as the corner.
  if (kmax < 0)
    reg_c = reg_param(lr); 
    rho_c = rho(lr); 
    eta_c = eta(lr);
  else
    index = find(lrho < x_corner & leta < y_corner);
    if ~isempty(index)
      [~,rpi] = min((lrho(index)-x_corner).^2 + (leta(index)-y_corner).^2);
      rpi = index(rpi);
    else
      [~,rpi] = min((lrho-x_corner).^2 + (leta-y_corner).^2);
    end
    reg_c = reg_param(rpi); 
    rho_c = rho(rpi); 
    eta_c = eta(rpi);
  end
  selection=0;

else
  error('Illegal method')
end
end