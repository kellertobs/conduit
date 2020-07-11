
if plot_op
    % prepare for plotting
    TX = {'Interpreter','Latex'}; FS = {'FontSize',16};
    TL = {'TickLabelInterpreter','Latex'}; TS = {'FontSize',12};
    UN = {'Units','Centimeters'};
    
    % set axis and border dimensions
    axh = 10.00; axw = 4*axh*R/L;
    ahs = 1.00; avs = 1.00;
    axb = 0.75; axt = 1.00;
    axl = 1.75; axr = 1.00;
    
    % initialize figures and axes
    fh1 = figure(1); clf; colormap(ocean);
    fh = axb + 1*axh + 0*avs + axt;
    fw = axl + 3*axw + 2*ahs + axr;
    set(fh1,UN{:},'Position',[5 5 fw fh]);
    set(fh1,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
    set(fh1,'Color','w','InvertHardcopy','off');
    set(fh1,'Resize','off');
    ax(11) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
    ax(12) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);
    ax(13) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+0*axh+0*avs axw axh]);

    fh2 = figure(2); clf; colormap(ocean);
    fh = axb + 1*axh + 0*avs + axt;
    fw = axl + 3*axw + 2*ahs + axr;
    set(fh2,UN{:},'Position',[10 10 fw fh]);
    set(fh2,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
    set(fh2,'Color','w','InvertHardcopy','off');
    set(fh2,'Resize','off');
    ax(21) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
    ax(22) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);
    ax(23) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+0*axh+0*avs axw axh]);
    
    fh3 = figure(3); clf; colormap(ocean);
    fh = axb + 1*axh + 0*avs + axt;
    fw = axl + 4*axw + 3*ahs + axr;
    set(fh3,UN{:},'Position',[15 15 fw fh]);
    set(fh3,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
    set(fh3,'Color','w','InvertHardcopy','off');
    set(fh3,'Resize','off');
    ax(31) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
    ax(32) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);
    ax(33) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+0*axh+0*avs axw axh]);
    ax(34) = axes(UN{:},'position',[axl+3*axw+3*ahs axb+0*axh+0*avs axw axh]);

    if plot_cv
        fh4 = figure(4); clf; colormap(ocean);
        fh = axb + 1*axh + 0*avs + axt;
        fw = axl + 4*axw + 3*ahs + axr;
        set(fh4,UN{:},'Position',[20 20 fw fh]);
        set(fh4,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
        set(fh4,'Color','w','InvertHardcopy','off');
        set(fh4,'Resize','off');
        ax(41) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
        ax(42) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);
        ax(43) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+0*axh+0*avs axw axh]);
        ax(44) = axes(UN{:},'position',[axl+3*axw+3*ahs axb+0*axh+0*avs axw axh]);
    end
    
    % plot velocity-pressure solution in Fig. 1
    figure(1);
    axes(ax(11));
    imagesc(r(2:end-1),z(2:end-1),-W(:      ,2:end-1)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$W$ [m/s]'],TX{:},FS{:});
    axes(ax(12));
    imagesc(r(2:end-1),z(2:end-1), U(2:end-1,:      )); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$U$ [m/s]'],TX{:},FS{:}); set(gca,'YTickLabel',[]);
    axes(ax(13));
    imagesc(r(2:end-1),z(2:end-1), P(2:end-1,2:end-1)./1e3); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$P$ [kPa]'],TX{:},FS{:}); set(gca,'YTickLabel',[]);
    drawnow;
    
    % plot phase-temperature solution in Fig. 2
    figure(2);
    axes(ax(21));
    imagesc(r(2:end-1),z(2:end-1), f(2:end-1,2:end-1) ); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\phi$ [vol]'],TX{:},FS{:});
    axes(ax(22));
    imagesc(r(2:end-1),z(2:end-1), c(2:end-1,2:end-1) ); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\chi$ [vol]'],TX{:},FS{:}); set(gca,'YTickLabel',[]);
    axes(ax(23));
    imagesc(r(2:end-1),z(2:end-1), T(2:end-1,2:end-1)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$T \ [^\circ$C]'],TX{:},FS{:}); set(gca,'YTickLabel',[]);
    drawnow;
    
    % plot density and rheology fields in Fig. 3
    figure(3);
    axes(ax(31));
    imagesc(r(2:end-1),z(2:end-1),      rho(2:end-1,2:end-1) ); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\rho$ [kg/m$^3$]'],TX{:},FS{:}); set(gca,'YTickLabel',[]);
    axes(ax(32));
    imagesc(r(2:end-1),z(2:end-1),log10(eta(2:end-1,2:end-1))); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\eta$ [log$_{10}$ Pas]'],TX{:},FS{:}); set(gca,'YTickLabel',[]);
    axes(ax(33));
    imagesc(r(2:end-1),z(2:end-1),log10(eII(2:end-1,2:end-1))); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\varepsilon_{II}$ [log$_{10}$ 1/s]'],TX{:},FS{:});
    axes(ax(34));
    imagesc(r(2:end-1),z(2:end-1),log10(tII(2:end-1,2:end-1))); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\tau_{II}$ [log$_{10}$ Pa]'],TX{:},FS{:}); set(gca,'YTickLabel',[]);
    drawnow;
    
    if plot_cv
        % plot residual fields in Fig. 4
        figure(4);
        axes(ax(41));
        imagesc(r(2:end-1),z(2:end-1),res_W.*dtW   ./(1e-16+norm(W(:),2)./N)); axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['res. $W$'],TX{:},FS{:});
        axes(ax(42));
        imagesc(r(2:end-1),z(2:end-1),res_U.*dtU   ./(1e-16+norm(U(:),2)./N)); axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['res. $U$'],TX{:},FS{:}); set(gca,'YTickLabel',[]);
        axes(ax(43));
        imagesc(r(2:end-1),z(2:end-1),res_P.*dtP   ./(1e-16+norm(P(:),2)./N)); axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['res. $P$'],TX{:},FS{:}); set(gca,'YTickLabel',[]);
        axes(ax(44));
        imagesc(r(2:end-1),z(2:end-1),res_f.*dt/10 ./(1e-16+norm(f(:),2)./N)); axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['res. $\phi$'],TX{:},FS{:}); set(gca,'YTickLabel',[]);
        drawnow;
    end
end

% save output to file
if save_op
    name = ['../out/',runID,'/',runID,'_vp_',num2str(floor(step/nop))];
    print(fh1,name,'-dpng','-r300','-opengl');
    name = ['../out/',runID,'/',runID,'_ft_',num2str(floor(step/nop))];
    print(fh2,name,'-dpng','-r300','-opengl');
    name = ['../out/',runID,'/',runID,'_mt_',num2str(floor(step/nop))];
    print(fh3,name,'-dpng','-r300','-opengl');
    
    name = ['../out/',runID,'/',runID,'_',num2str(floor(step/nop))];
    save(name,'U','W','P','f','T','c','dfdt','dTdt','rho','CL','eta','Div_V','err','ezz','erz','trr','tzz','trz','eII','tII','dt','time','step','fvol0');
    name = ['../out/',runID,'/',runID,'_cont'];
    save(name,'U','W','P','f','T','c','dfdt','dTdt','rho','CL','eta','Div_V','err','ezz','erz','trr','tzz','trz','eII','tII','dt','time','step','fvol0');
    
    if step == 1
        logfile = ['../out/',runID,'/',runID,'.log'];
        if exist(logfile,'file'); delete(logfile); end
        diary(logfile)
    end
end
    