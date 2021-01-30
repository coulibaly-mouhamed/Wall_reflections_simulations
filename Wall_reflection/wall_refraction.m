function[theta_r] =wall_refraction(theta)

% Input: theta -- incident angle

    %% Setup parameters
    mem = 0.9;
    Gam = mem*4.20;
    Nx = 128; Ny = Nx; 
    Lx = 32; Ly = Lx; dt_desired = min(Lx/Nx,Ly/Ny)/8;
    plotoption=inf;
    p = problem_setup_refraction(Nx,Ny,Lx,Ly,Gam,dt_desired);

    %% Initial position and velocity
    theta_rad = theta*pi/180;
    vfreespace = 0.0555;
    dist = 15;
    p.xi = 10-dist; p.yi = mod(-dist*tan(theta_rad)+p.Ly/2,p.Ly)-p.Ly/2; 
    p.ui=0.1*cos(theta_rad); p.vi=0.1*sin(theta_rad);
    p.nimpacts = min(round(2*(dist./cos(theta_rad)./vfreespace)),2400);
    p.nimpacts = 900;

    %% Bottom profile
    p.h   = p.h0_shallow.*(p.xx>10)+p.h0_deep.*(p.xx<=10);
    p.d   = p.d0_shallow.*(p.xx>10)+p.d0_deep.*(p.xx<=10);
    [x_data,y_data,t_data,eta_data] = trajectory(p,plotoption);
    n = length(x_data);
    u_m =[0,10];
    u_d = [x_data(n)-x_data(n-1), y_data(n)-y_data(n-1)];
    theta_r =acos(dot(u_m, u_d)/(norm(u_m)*norm(u_d))); %peut etre faire une moyenne serait plus interessant...
    %plot(x_data,y_data);
    p = gather(p);

    %save(['reflection_wavefield_thetaI_',num2str(theta,'%.2f'),'_r_',num2str(p.drop_radius),...
     %     '_Gam_',num2str(p.Gam),'_N_',num2str(p.Nx),'_L_',num2str(Lx),'.mat'],...
      %      'x_data','y_data','t_data','p','eta_data');
    
end
