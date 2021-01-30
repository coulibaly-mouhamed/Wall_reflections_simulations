function plot_solution(eta,xi,yi,t,p,cplot)
    figure(1)
    warning off;
    if p.useGPU ==1
        eta = gather(eta);
    end
    
    if nargin < 6
        cplot=[-0.001,0.002];
    end
    eta_interp = interp2(p.xx,p.yy,eta,p.xxx,p.yyy,'spline',0);
    h = pcolor(p.xxx,p.yyy,eta_interp);set(h,'edgecolor','none','FaceAlpha',0.95); grid off; 
    colorbar; caxis(cplot)
    title(['t=',num2str(t)]); hold on;
    switch p.num_drops
        case 1
        plot(xi,yi,'k.','MarkerSize',10)
        v=[p.d0_shallow-0.0001,p.d0_shallow+0.0001];
        contour(p.xx,p.yy,p.d,v,'LineWidth',2,'LineColor','k');
        case 2
        plot(xi(1),yi(1),'g.','MarkerSize',10)
        plot(xi(2),yi(2),'r.','MarkerSize',10)
        v=[p.d0_shallow-0.0001,p.d0_shallow+0.0001];
        contour(p.xx,p.yy,p.d,v,'LineWidth',2,'LineColor','k');
    end
    axis square;
    drawnow;  hold off; 

    
if p.plotPS == 1
    figure(2)
    clf
    PS = abs(fft2(eta));
    [C,h2] = contourf(imag(p.Kx)/(2*pi),imag(p.Ky)/(2*pi),PS/max(max(PS)));
    set(h2,'LineStyle','none');
    colorbar;
    axis square;
    drawnow;
end
    
    