function plot_wavefield_crossection(eta,xi,yi,t,p)
    figure(1); %clf;    
    plotwindow = 5;
    xxx = linspace(xi-plotwindow,xi+plotwindow,500); yyy = linspace(yi-plotwindow,yi+plotwindow,500);
    eta_interp = interp2(p.xx,p.yy,eta,xxx,yi,'spline',0);
    plot(xxx,eta_interp); 
    title(['t=',num2str(t)]); hold on;

    xlim([xi-plotwindow,xi+plotwindow]); 
    xlabel('x'); ylabel('\eta');
    %axis square;
    grid on; drawnow;  
    
    figure(2); clf;    
    eta_interp = interp2(p.xx,p.yy,eta,xi,yyy,'spline',0);
    plot(yyy,eta_interp);
    title(['t=',num2str(t)]); hold on;
    xlim([yi-plotwindow,yi+plotwindow]); 
    xlabel('y'); ylabel('\eta');
    %axis square;
    grid on;drawnow; 
    
   