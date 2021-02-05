function [theta_r,y_a,y_b] =wall_reflection(theta)

% Input: theta -- incident angle

    %% Setup parameters
    mem = 0.9;
    Gam = mem*4.20;
    Nx = 128; Ny = Nx; 
    Lx = 32; Ly = Lx; dt_desired = min(Lx/Nx,Ly/Ny)/8;
    plotoption=40;
    p = problem_setup_reflection(Nx,Ny,Lx,Ly,Gam,dt_desired);

    %% Initial position and velocity
    theta_rad = theta*pi/180;
    vfreespace = 0.0555;
    dist =2;
    p.xi = 2-dist; p.yi = mod(-dist*tan(theta_rad)+p.Ly/2,p.Ly)-p.Ly/2; 
    p.ui=0.1*cos(theta_rad); p.vi=0.1*sin(theta_rad);
    %p.nimpacts = min(round(2*(dist./cos(theta_rad)./vfreespace)),2400);
    p.nimpacts=700;

    %% Bottom profile
    p.h   = p.h0_shallow.*(p.xx>10)+p.h0_deep.*(p.xx<=10);
    p.d   = p.d0_shallow.*(p.xx>10)+p.d0_deep.*(p.xx<=10);
    disp(num2str(p.nimpacts));
    p.makeMovie=0;
    [x_data,y_data,t_data,eta_data] = trajectory(p,plotoption);
    n = length(x_data);
    u_m =[1,0];
    u_d = [x_data(n)-x_data(n-1), y_data(n)-y_data(n-1)];
    %calcul des poits A et B
    %Pour A
    a_a= y_data(2)-y_data(1)/x_data(2)-x_data(1);
    b_a = y_data(2)-a_a*x_data(2);
    y_a = a_a*10+b_a;
    %Pour B
    a_b= y_data(n)-y_data(n-1)/x_data(n)-x_data(n-1);
    b_b = y_data(n)-a_b*x_data(n);
    y_b = a_b*10+b_b;
    theta_r =acos(dot(u_m, u_d)/(norm(u_m)*norm(u_d))); %peut etre faire une moyenne serait plus interessant...
    p = gather(p);

    %save(['reflection_wavefield_thetaI_',num2str(theta,'%.2f'),'_r_',num2str(p.drop_radius),...
     %     '_Gam_',num2str(p.Gam),'_N_',num2str(p.Nx),'_L_',num2str(Lx),'.mat'],...
      %      'x_data','y_data','t_data','p','eta_data');
    
end
