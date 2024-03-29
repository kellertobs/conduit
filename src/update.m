% update solution-dependent parameter and auxiliary variable fields
eII(2:end-1,2:end-1) = 1e-16 + (0.5.*(err(2:end-1,2:end-1).^2 + ezz(2:end-1,2:end-1).^2 ...  % get strain rate magnitude
    + 2.*(erz(1:end-1,1:end-1).^2.*erz(2:end,1:end-1).^2.*erz(1:end-1,2:end).^2.*erz(2:end,2:end).^2).^0.25)).^0.5;  
eII(:,[1 end]) = eII(:,[end-1 2]);
eII([1 end],:) = eII([2 end-1],:);

tII(2:end-1,2:end-1) = 1e-16 + (0.5.*(trr(2:end-1,2:end-1).^2 + tzz(2:end-1,2:end-1).^2 ...  % get stress magnitude
    + 2.*(trz(1:end-1,1:end-1).^2.*trz(2:end,1:end-1).^2.*trz(1:end-1,2:end).^2.*trz(2:end,2:end).^2).^0.25)).^0.5;  
tII(:,[1 end]) = tII(:,[end-1 2]);
tII([1 end],:) = tII([2 end-1],:);

eta   = eta0 .* max(1e-6,1-2*f).^-A .* max(1e-6,1-2*c).^-B;                % get bubble-crystal-dep. magma viscosity
eta   = 1./(1./(eta + eta0./sqrt(etactr)) + 1./(eta0.*sqrt(etactr)));
for d = 1:ceil(delta)
    dd  = delta/ceil(delta);
    eta = eta + dd.*(diff(eta([end-2,1:end,3],:),2,1)./8 + diff(eta(:,[1,1:end,end]),2,2)./8);
end
etac  = (eta(1:end-1,1:end-1)+eta(2:end,1:end-1) ...                       % get viscosity in cell corners
       + eta(1:end-1,2:end  )+eta(2:end,2:end  ))./4;
    
rho   = (1-f-c).*rhom + c.*rhoc + f.*rhof;                                 % get bubble-crystal-dep. magma density
rhoBF = (rho(1:end-1,2:end-1)+rho(2:end,2:end-1))./2-mean(rho(:));         % get relative density for bouancy force term

C     = (1-f-c).*Cm + c.*Cc + f.*Cf;                                       % get bubble-crystal-dep. magma heat capacity
CL    = C + (1-f).*LH./(Tliq-Tsol) .* double(T>Tsol & T<Tliq);             % get heat capacity adjusted for latent heat

dtW   = (h/2)^2./((1+zeta).*max(eta(1:end-1,:),eta(2:end,:)));             % iterative step size
dtU   = (h/2)^2./((1+zeta).*max(eta(:,1:end-1),eta(:,2:end)));             % iterative step size
dtP   = (2/3+zeta).*eta;                                                 % iterative step size
dtT   = min((h/2)^2./kT.*rho(:).*CL(:));                                   % diffusive time step size
dta   = min(CFL*min(h/max(abs([U(:);W(:)]+1e-16))/2));                     % advective time step size
dt    = min(2.*dt,min(dta,dtT));                                           % physical time step size