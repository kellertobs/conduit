% get residual norm
resnorm = norm(res_U(:).*dtU(:),2)./(1e-16+norm(U(:),2)) ...
        + norm(res_W(:).*dtW(:),2)./(1e-16+norm(W(:),2)) ...
        + norm(res_P(:).*dtP(:),2)./(1e-16+norm(P(:),2)) ...
        + norm(res_f(:).*dt/10 ,2)./(1e-16+norm(f(:),2)) ...
        + norm(res_T(:).*dt/10 ,2)./(1e-16+norm(T(:),2));

if it==0 || resnorm>resnorm0; resnorm0 = resnorm; end  % set reference residual

% get total mass error
fvol = sum(f(:));
if step == 1 && it<=nup; fvol0 = fvol; end  % set reference residual
ferr  = abs(fvol-fvol0)/fvol0;

% report iterations
if     it >=  0  && it <  10
    fprintf(1,'    ---  it =      %d;   abs res = %4.4e;   rel res = %4.4e;   err = %4.4e \n',it,resnorm,resnorm/resnorm0,ferr);
elseif it >= 10  && it < 100
    fprintf(1,'    ---  it =     %d;   abs res = %4.4e;   rel res = %4.4e;   err = %4.4e \n',it,resnorm,resnorm/resnorm0,ferr);
elseif it >= 100 && it < 1000
    fprintf(1,'    ---  it =    %d;   abs res = %4.4e;   rel res = %4.4e;   err = %4.4e \n',it,resnorm,resnorm/resnorm0,ferr);
elseif it >= 1000 && it < 10000
    fprintf(1,'    ---  it =   %d;   abs res = %4.4e;   rel res = %4.4e;   err = %4.4e \n',it,resnorm,resnorm/resnorm0,ferr);
elseif it >= 10000
    fprintf(1,'    ---  it =  %d;   abs res = %4.4e;   rel res = %4.4e;   err = %4.4e \n',it,resnorm,resnorm/resnorm0,ferr);
end 

% plot convergence of outer iterations
if plot_cv
    figure(100); if it==0; clf; else; hold on; end
    plot(it,log10(resnorm),'r.','MarkerSize',15); box on; axis tight;
    drawnow;
end

% output