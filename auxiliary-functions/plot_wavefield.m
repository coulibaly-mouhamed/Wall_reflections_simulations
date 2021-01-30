function plot_wavefield(eta,xi,yi,t,p,cplot)
    figure(1)
    
    %% Convert to dimensional units
    
    eta = (0.005)*10^(6)*eta; % Eta in micrometers
    p.xx = p.xF*1000 * p.xx; p.yy =p.xF*1000*p.yy;  % Space in milimeters
    xi = p.xF*1000* xi; yi =p.xF*1000*yi;
    wall_position = p.xF*1000*10;
    
    if nargin < 6
        cplot=[min(min(eta)),max(max(eta))];
    end
    
    %subplot(1,2,1)
    plotwindow = 25;
    %xxx = linspace(xi-plotwindow,xi+plotwindow,500); yyy = linspace(yi-plotwindow,yi+plotwindow,500);
    %xxx = linspace(-1,15,500); yyy = linspace(-8,8,500);
    xxx = linspace(-150,150,500); yyy = linspace(-150,150,500);
    
    [XX,YY] = meshgrid(xxx,yyy); 
    %eta = real(ifft2(fft2(eta).*(abs(p.K2)<1.2*(2*pi)^2))); % Possible filter for high frequencies
    eta_interp = interp2(p.xx,p.yy,eta,XX,YY,'spline',0); % Prefactor converts to micrometers
    colormap parula;
    %colormap autumn;
    h = pcolor(xxx,yyy,eta_interp);shading interp;
    %set(h,'edgecolor','none','FaceAlpha',0.95); 
    %view(45,85);
    grid off; 
    cb = colorbar; 
    ylabel(cb,'$\eta \ (\mu m)$','interpreter','latex');
    caxis(cplot)
    light
    %             camlight
    lighting phong
    material([0.85 0.4 0.1 1.5])

    title(['t=',num2str(t*0.025,'%.3f'),' s']); hold on;
    plot3(xi,yi,interp2(p.xx,p.yy,eta,xi,yi,'spline'),'k.','MarkerSize',20);
    plot([wall_position,wall_position],[-100,100],'k','LineWidth',3);
    xlim([xxx(1),xxx(end)]); ylim([yyy(1),yyy(end)]);
    xlabel('x (mm)'); ylabel('y (mm)');
    axis square;
    drawnow;  hold off; 
    
%     subplot(1,2,2)
%     plot(xi,yi,'k.','MarkerSize',20);
%     shading interp; 
%     plot([10,10],[-p.Ly/2,p.Ly/2],'k');
%     title(['t=',num2str(t)]); hold on;
%     xlim([-p.Lx/2,p.Lx/2]); ylim([-p.Ly/2,p.Ly/2]); 
%     xlabel('x'); ylabel('y');
%     axis square;
    
    