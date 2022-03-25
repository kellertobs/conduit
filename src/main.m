% get coordinate arrays
r     = -h/2:h:R+h/2;
z     = -h/2:h:L+h/2;
[rr,zz] = meshgrid(r,z);
rc    = (r(1:end-1)+r(2:end))./2;
zc    = (z(1:end-1)+z(2:end))./2;

% get smoothed initialisation field
rng(15);
a = randn(N,(N-2)*R/L+2);
b = zeros(N,(N-2)*R/L+2); b(rr<R/2 & abs(zz-L/2)<L/(2*SlugNo)) = SlugNo;
for i = 1:round(smth)
    a(2:end-1,2:end-1) = a(2:end-1,2:end-1) + diff(a(:,2:end-1),2,1)./8 + diff(a(2:end-1,:),2,2)./8;
    a([1 end],:) = a([end-1 2],:);
    a(:,[1 end]) = a(:,[2 end-1]);
    b(2:end-1,2:end-1) = b(2:end-1,2:end-1) + diff(b(:,2:end-1),2,1)./8 + diff(b(2:end-1,:),2,2)./8;
    b([1 end],:) = b([end-1 2],:);
    b(:,[1 end]) = b(:,[2 end-1]);
end
a = a./max(abs(a(:)));
b = b./max(abs(b(:))).*SlugNo;

% set initial solution fields
f      =  max(1e-3,min(1-1e-3,(1+a).*f1 + b.*f2));  fo = f;  fi = f;  res_f = 0.*f;  fvol0 = sum(f(:));
T      =  T0.*ones(size(f));  To = T;  Ti = T;  res_T = 0.*T;
c      =  max(0,min(1-f,(1-f).*(T-Tliq)./(Tsol-Tliq)));
U      =  zeros(size((rr(:,1:end-1)+rr(:,2:end))));  Ui = U;  res_U = 0.*U;
W      =  zeros(size((rr(1:end-1,:)+rr(2:end,:))));  Wi = W;  res_W = 0.*W;
P      =  0.*f;  Pi = P;  res_P = 0.*P;  meanQ = 0;

% initialise auxiliary variables
dfdt   =  0.*f;  lapl_f = 0.*f;  
dTdt   =  0.*T;  lapl_T = 0.*T; 
advn   =  0.*f;
eIIref =  1e-6;  
Div_V  =  0.*P;
err    =  0.*P;  ezz = 0.*P;  erz = zeros(size(f)-1);  eII = 0.*P;  
trr    =  0.*P;  tzz = 0.*P;  trz = zeros(size(f)-1);  tII = 0.*P; 

% initialise timing and iterative parameters
step   =  0;
time   =  0;
dt     =  dt/2;
it     =  0;

% overwrite fields from file if restarting run
if     restart < 0  % restart from last continuation frame
    name = ['../out/',runID,'/',runID,'_cont'];
    load(name,'U','W','P','f','T','c','dfdt','dTdt','rho','CL','eta','Div_V','err','ezz','erz','trr','tzz','trz','eII','tII','dt','time','step','fvol0');
    name = ['../out/',runID,'/',runID,'_par'];
    load(name);
elseif restart > 0  % restart from specified continuation frame
    name = ['../out/',runID,'/',runID,'_',num2str(restart)];
    load(name,'U','W','P','f','T','c','dfdt','dTdt','rho','CL','eta','Div_V','err','ezz','erz','trr','tzz','trz','eII','tII','dt','time','step','fvol0');
    name = ['../out/',runID,'/',runID,'_par'];
    load(name);
end

load ocean;  % load custom colormap

% initialise material properties and plot initial condition
update;  output;

% physical time stepping loop
while time <= tend && step <= M
        
    fprintf(1,'\n\n\n*****  step %d;  dt = %4.4e;  time = %4.4e;\n\n',step,dt,time);
    tic;
    
    % store previous solution
    fo      = f;
    co      = c;
    To      = T;
    dfdto   = dfdt;
    dTdto   = dTdt;
    
    % reset residuals and iteration count
    resnorm  = 1e3;
    resnorm0 = resnorm;
    it       = 0;
    
    % non-linear iteration loop
    startup = 2*double(step<=0) + double(step>0);
    while resnorm/resnorm0 >= rtol^startup && resnorm >= atol/startup && it <= maxit*startup || it <= nup
                
        % store previous iterative solutions   
        Uii = Ui;  Ui = U;
        Wii = Wi;  Wi = W;
        Pii = Pi;  Pi = P;
        fii = fi;  fi = f;
        Tii = Ti;  Ti = T;

        % update fields
        if ~mod(it,nup); update; end

        if ~mod(it,nup)
            % update bubble fraction
            a = f; advection; VGrd_f = advn; %#ok<NASGU>                   % get upwind advection term
            
            qfz   = - kf.*(f(1:end-1,:)+f(2:end,:))/2 .* ddz(f,h);         % z-volume diffusion flux
            qfr   = - kf.*(f(:,1:end-1)+f(:,2:end))/2 .* ddr(f,h);         % r-volume diffusion flux
            lapl_f(2:end-1,2:end-1) = - ddz(qfz(:,2:end-1),h) ...          % laplacian term
                                      - 1./r(2:end-1) .* ddr(rc.*qfr(2:end-1,:),h);  
            
            dfdt  = lapl_f - VGrd_f;                                       % total rate of change
            
            res_f = (f-fo)./dt - theta.*dfdt - (1-theta).*dfdto;           % residual bubble evolution
         
            res_f([1,end],:) = res_f([end-1 2],:);                         % periodic top/bot boundaries
            res_f(:,[1 end]) = res_f(:,[2 end-1]);                         % zero-flux side boundaries
            
            res_f = min(res_f,fi./alpha./dt*4);                            % prevent f < 0

            f = fi - alpha.*res_f.*dt/10 + beta.*(fi-fii);                  % update bubble solution
            
            % update temperature
            a = T; advection; VGrd_T = advn;                               % get upwind advection term
            
            qTz    = - (kT(1:end-1,:)+kT(2:end,:))./2 .* ddz(T,h);                                     % z-heat diffusion flux
            qTr    = - (kT(:,1:end-1)+kT(:,2:end))./2 .* ddr(T,h);                                     % r-heat diffusion flux
            lapl_T(2:end-1,2:end-1) = (- ddz(qTz(:,2:end-1),h) ...         % laplacian term
                                       - 1./r(2:end-1) .* ddr(rc.*qTr(2:end-1,:),h)) ...
                                       ./rho(2:end-1,2:end-1)./CL(2:end-1,2:end-1);  
                                  
            cool   =  (T-100)./tau_c.*exp((rr-R)/4/h);                     % get constant-flux wall cooling, inflow heating
            heat   = -(T-T0 )./dt/4 .*exp((zz-L)/4/h) .* (f>f1 & (W([end-1,1:end],:)+W([1:end,2],:))./2<=0);
            
            dTdt   = lapl_T  - VGrd_T - cool + heat;                       % total rate of change
            res_T  = (T-To)./dt - theta.*dTdt - (1-theta).*dTdto;          % residual temperature evolution         
            
            res_T([1,end],:) = res_T([end-1 2],:);                         % periodic top/bot boundaries
            res_T(:,[1 end]) = res_T(:,[2 end-1]);                         % zero-flux side boundaries
            
            T = Ti - alpha.*res_T.*dt/10 + beta.*(Ti-Tii);                 % update temperature solution

            c = max(0,min(1-f,(1-f).*(T-Tliq)./(Tsol-Tliq)));              % update crystallinity
            
            if step==1; dfdto = dfdt;  dTdto = dTdt; end
        end
        
        % update z-velocity
        Div_tz = ddz(tzz(:,2:end-1),h) + ddr(trz,h);                       % get z-stress divergence
        cyl_tz = (trz(:,1:end-1)+trz(:,2:end))./2./r(:,2:end-1);           % cylindrical z-stress term
        
        meanW = sum(sum(r(2:end-1).*W(:,2:end-1)))./sum(sum(r(2:end-1).*ones(size(W(:,2:end-1)))));
        
        res_W(:,2:end-1) = - Div_tz - cyl_tz  ...                          % residual z-momentum conservation
                           + ddz(P(:,2:end-1),h) - rhoBF.*g0  ...
                           + meanW./dtW(:,2:end-1);                         

        res_W([1 end],:) = [sum(res_W([1 end],:),1)./2; ...
                            sum(res_W([1 end],:),1)./2];                   % periodic boundaries
        res_W(:,1  ) =  res_W(:,2    );                                    % free slip inner side boundary
        res_W(:,end) = -res_W(:,end-1);                                    % no slip outer side boundary
                
        W = Wi - alpha.*res_W.*dtW + beta.*(Wi-Wii);                       % update z-velocity solution

        % update r-velocity        
        Div_tr = ddr(trr(2:end-1,:),h) + ddz(trz,h);                       % get r-stress divergence
        cyl_tr = (etac(1:end-1,:)+etac(2:end,:))./2./(rc+1e-3) .* ...      % cylindrical z-stress term
                 (ddr((U(2:end-1,[1,1:end])+U(2:end-1,[1:end,end]))./2,h) ...
                 - U(2:end-1,:)./(rc+1e-3));                                      
        
        meanU = sum(sum(rc.*U(2:end-1,:)))./sum(sum(rc.*ones(size(U(2:end-1,:)))));

        res_U(2:end-1,:) = - Div_tr - cyl_tr ... 
                           + ddr(P(2:end-1,:),h);                          % residual r-momentum conservation
        
        res_U([1 end],:) = res_U([end-1 2],:);                             % periodic top/bot boundaries
        res_U(:,[1 end]) = 0;                                              % no flow across side boundaries


        U = Ui - alpha.*res_U.*dtU + beta.*(Ui-Uii);                       % update r-velocity solution
        
        % update velocity divergence
        Div_V(2:end-1,2:end-1) = ddz(W(:,2:end-1),h) ...                   % get velocity divergence
                               + 1./r(2:end-1) .* ddr(rc.*U(2:end-1,:),h);                      
        Div_V([1 end],:) = Div_V([end-1 2],:);                             % periodic top/bot boundaries
        Div_V(:,[1 end]) = Div_V(:,[2 end-1]);                             % continuous side boundaries
        
        % update strain rates
        err(:,2:end-1)   = ddr(U,h) - Div_V(:,2:end-1)./3;               % get r-normal strain rate
        err([1 end],:)   = err([end-1 2],:);                               % periodic top/bot boundaries
        err(:,[1 end])   = err(:,[2 end-1]);                               % continuous side boundaries
        ezz(2:end-1,:)   = ddz(W,h) - Div_V(2:end-1,:)./3;               % get z-normal strain rate
        ezz([1 end],:)   = ezz([end-1 2],:);                               % periodic top/bot boundaries
        ezz(:,[1 end])   = ezz(:,[2 end-1]);                               % continuous side boundaries
        erz              = 1/2.*(ddz(U,h) + ddr(W,h));                     % get shear strain rate

        % update stresses
        trr = eta .* err;                                                  % r-normal stress
        tzz = eta .* ezz;                                                  % z-normal stress
        trz = etac.* erz;                                                  % rz-shear stress  
        
        % update dynamic pressure        
        res_P = Div_V;                                                     % residual mass conservation
        meanP = sum(sum(r(2:end-1).*P(2:end-1,2:end-1)))./sum(sum(r(2:end-1).*ones(size(P(2:end-1,2:end-1)))));

        res_P([1 end],:) = res_P([end-1 2],:);                             % periodic top/bot boundaries
        res_P(:,[1 end]) = res_P(:,[2 end-1]);                             % continuous side boundaries
        
        P = Pi - alpha.*res_P.*dtP + beta.*(Pi-Pii);                       % update pressure solution
        
        % check and report convergence every nup iterations
        if ~mod(it,nup); report; end
        
        it = it+1;
    end
    
    % print diagnostics
    fprintf(1,'\n         time to solution = %4.4f sec\n\n',toc);
    
    fprintf(1,'         min f   =  %1.4f;    mean f   = %1.4f;    max f   = %1.4f;   [vol]\n'  ,min(min(f(2:end-1,2:end-1)  )),mean(mean(f(2:end-1,2:end-1)  ))   ,max(max(f(2:end-1,2:end-1)  )));
    fprintf(1,'         min c   =  %1.4f;    mean c   = %1.4f;    max c   = %1.4f;   [vol]\n'  ,min(c(:)  ),mean(c(:)  )   ,max(c(:)  ));
    fprintf(1,'         min T   =  %4.1f;    mean T   = %4.1f;    max T   = %4.1f;   [degC]\n\n' ,min(T(:)  ),mean(T(:)  )   ,max(T(:)  ));

    fprintf(1,'         min rho =  %4.1f;    mean rho = %4.1f;    max rho = %4.1f;   [kg/m3]\n'  ,min(rho(:)),mean(rho(:))   ,max(rho(:)));
    fprintf(1,'         min eta =  %1.2e;  mean eta = %1.2e;  max eta = %1.2e; [Pas]\n\n',min(eta(:)),geomean(eta(:)),max(eta(:)));

    fprintf(1,'         min U   = %1.4f;    mean U   = %1.4f;    max U   = %1.4f;   [m/s]\n'  ,min(U(:)  ),mean(U(:)  ),max(U(:)  ));
    fprintf(1,'         min W   = %1.4f;    mean W   = %1.4f;    max W   = %1.4f;   [m/s]\n'  ,min(-W(:) ),mean(-W(:) ),max(-W(:) ));
    fprintf(1,'         min P   = %2.4f;   mean P   = %2.4f;    max P   = %2.4f;  [kPa]\n\n',min(P(:)./1e3),mean(P(:)./1e3),max(P(:)./1e3));

    % plot results
    if ~mod(step,nop); output; end
    
    % increment time/step
    time = time+dt;
    step = step+1;
    
end

diary off
