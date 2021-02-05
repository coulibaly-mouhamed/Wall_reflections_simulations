function [x_data,y_data,t_data,eta_data,dist] = trajectory(p,plotgap)

% Set initial condition 
t=0;
phi     = p.phi0;
eta     = p.eta0; 
xi = p.xi; yi=p.yi; ui=p.ui; vi=p.vi; 

dist =[];
% Create vector for position
x_data = zeros(p.nimpacts,p.num_drops); y_data = x_data; t_data=zeros(p.nimpacts,1);

% Create vector for wavefield
if p.store_wavefield == 0; 
    eta_data = nan;
else
    eta_data = zeros(p.Nx,p.Ny,length(p.store_wavefield)); 
end

% % Move wave field to GPU
if p.useGPU == 1
    eta = gpuArray(eta);
    phi = gpuArray(phi);
    x_data = gpuArray(x_data);
    y_data = gpuArray(y_data);
    t_data = gpuArray(t_data);
end

% Fourier transform (2d) of initial condition
phi_hat = fft2(phi);               
eta_hat = fft2(eta);
ii = 1;
switch p.makeMovie
    case 1
        vidfile = VideoWriter('wall_reflection.mp4','MPEG-4');
        open(vidfile);
end
for n=1:p.nimpacts
    
    eta = real(ifft2(eta_hat));   
    eta_max=max(max(abs(eta)));
    if eta_max > 1
        eta_max
        disp('Probably above threshold. Stopping simulation...');
        plot_solution(eta,xi,yi,t,p);
        break
    end   
    
    % Store position (and possibly wavefield)
    x_data(n,:) = xi; y_data(n,:) = yi; t_data(n) = t;
   
    if ismember(n,p.store_wavefield) == 1;
        eta_data(:,:,ii) = gather(eta); 
        ii = ii+1;
    end
    
    % Plot
    switch p.makeMovie
        case 1
            if mod(n-1,plotgap)==0
                plot_solution(eta,xi,yi,t,p);
                F = getframe(gcf);
                title('Reflection');
                writeVideo(vidfile, F);
            end
        case 0
            if mod(n-1,plotgap)==0
                plot_solution(eta,xi,yi,t,p);
            end
    end

     % Drop impact
    [ui, vi, phi_hat] = drop_impact(xi,yi, ui, vi, phi_hat, eta_hat, p);

    % Evolve drops between impacts
    [xi, yi, ui, vi] = evolve_drops(xi, yi, ui, vi, p);
    % Add noise (if included)
    ui = ui + p.sig_noise.*randn; vi = vi + p.sig_noise.*randn;

    
    % Evolve wave between impacts
    [phi_hat, eta_hat]           = evolve_wave(phi_hat, eta_hat, t, p);
    switch p.num_drops
        case 2
            dist= [dist sqrt((xi(1)-xi(2)).^2+(yi(1)-yi(2)).^2)];
    end
    
        % Diplay results
    %if mod(n-1,100)==0
     %   ui_av = xi - x_data(n,:); vi_av = yi - y_data(n,:);
      %  disp(['impact=',num2str(n),'  xi=',num2str(xi),'  yi=',num2str(yi),...
       % '  ui=',num2str(ui_av),'  vi=',num2str(vi_av),...
        %'  |v|=',num2str(sqrt((ui_av).^2+(vi_av).^2)),'  eta_max=',num2str(eta_max),...
        %'  |eta_boundary|=',num2str(max(abs(eta(1,:))))]);
             %disp(['impact=',num2str(n)]);
    %end
    
    t = t+p.impact_interval;

end
if p.makeMovie==1
    close(vidfile)
end
p.dist =dist;
if p.useGPU ==1
    x_data = gather(x_data); y_data = gather(y_data); 
    eta_data = gather(eta_data);
    t_data = gather(t_data);
end

end





