% read out advected field at stencil positions
acc = a(2:end-1,2:end-1);
ajp = a(3:end  ,2:end-1);  ajpp = a([4:end,3      ],2:end-1);
ajm = a(1:end-2,2:end-1);  ajmm = a([end-2,1:end-3],2:end-1);
aip = a(2:end-1,3:end  );  aipp = a(2:end-1,[4:end,end    ]);
aim = a(2:end-1,1:end-2);  aimm = a(2:end-1,[1    ,1:end-3]);
        
% switch to chosen advection scheme
switch ADVN
    
    case 'UPW2'  % second-order upwind scheme
        
        % read out velocities at stencil positions
        wm   =  min((W(2:end  ,2:end-1)+W(1:end-1,2:end-1))./2,0);
        wp   =  max((W(2:end  ,2:end-1)+W(1:end-1,2:end-1))./2,0);
        um   =  min((U(2:end-1,2:end  )+U(2:end-1,1:end-1))./2,0);
        up   =  max((U(2:end-1,2:end  )+U(2:end-1,1:end-1))./2,0);
        
        % get upwind gradients of advected field
        da_dzm =   ( 3.*acc-4.*ajm+ajmm)./2./h;
        da_dzp =   (-3.*acc+4.*ajp-ajpp)./2./h;
        da_dxm =   ( 3.*acc-4.*aim+aimm)./2./h;
        da_dxp =   (-3.*acc+4.*aip-aipp)./2./h;
        
        % get advection term (v . Grad a)
        advn(2:end-1,2:end-1)  =  um .* da_dxp + up .* da_dxm + wm .* da_dzp + wp .* da_dzm;

    case 'UPW3'  % third-order upwind scheme
        
        % read out velocities at stencil positions
        wm   =  min((W(2:end  ,2:end-1)+W(1:end-1,2:end-1))./2,0);
        wp   =  max((W(2:end  ,2:end-1)+W(1:end-1,2:end-1))./2,0);
        um   =  min((U(2:end-1,2:end  )+U(2:end-1,1:end-1))./2,0);
        up   =  max((U(2:end-1,2:end  )+U(2:end-1,1:end-1))./2,0);
        
        % get upwind gradients of advected field
        da_dzm =   ( 2.*ajp+3.*acc-6.*ajm+ajmm)./6./h;
        da_dzp =   (-2.*ajm-3.*acc+6.*ajp-ajpp)./6./h;
        da_dxm =   ( 2.*aip+3.*acc-6.*aim+aimm)./6./h;
        da_dxp =   (-2.*aim-3.*acc+6.*aip-aipp)./6./h;
        
        % get advection term (v . Grad a)
        advn(2:end-1,2:end-1)  =  um .* da_dxp + up .* da_dxm + wm .* da_dzp + wp .* da_dzm;
        
    case 'FRM'  % Fromm upwind flux-conservative scheme
        
        % read out velocities at stencil positions
        wp = W(2:end  ,2:end-1);
        wm = W(1:end-1,2:end-1);
        up = U(2:end-1,2:end  );
        um = U(2:end-1,1:end-1);
        
        % get advection term (Div . a v - a Div . v)
        advn(2:end-1,2:end-1)  = ((up .*(-aipp + 5.*(aip+acc)-aim )./8 - abs(up).*(-aipp + 3.*(aip-acc)+aim )./8) ...
                               -  (um .*(-aip  + 5.*(acc+aim)-aimm)./8 - abs(um).*(-aip  + 3.*(acc-aim)+aimm)./8))./h ...
                               + ((wp .*(-ajpp + 5.*(ajp+acc)-ajm )./8 - abs(wp).*(-ajpp + 3.*(ajp-acc)+ajm )./8) ...
                               -  (wm .*(-ajp  + 5.*(acc+ajm)-ajmm)./8 - abs(wm).*(-ajp  + 3.*(acc-ajm)+ajmm)./8))./h ...
                               -   acc .* ((wp-wm)./h + (up-um)./h);

end