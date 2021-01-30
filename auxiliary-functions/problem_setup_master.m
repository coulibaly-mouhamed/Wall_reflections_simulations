function p = problem_setup_master(Nx, Ny, Lx, Ly, Gam, dt_desired)
% Sets most of the parameters for the problem
% Input: 
%   Nx          -------- Number of points in x
%   Ny          -------- Number of points in x
%   Lx          -------- Size of domain in x
%   Ly          -------- Size of domain in y
%   Gam         -------- Amplitude of the shaking
%   dt_desired  -------- Desired time step


%% Problem parameters in SI units (typically do not change for a given experiment)
nu        = 2*10^(-5);        % Kinematic viscosity (m2/s)
nu        = 0.8025*nu;        % Effective viscosity
rho       = 950;              % Density (kg/m3);
sig       = 0.0206;           % Surface tension (N/m);
h_deep    = 0.00609;            % depth  of deep region (m);
h_shallow = 0.00042;           % depth of shallow region (m)
omega0     = 80*2*pi;         % Angular frequency in 1/s
g0         = 9.8;             % Gravity in m/s
%Gam       = 4.2;  

%% Dispersion relation for inviscid surface waves
dispEuler = @(k,h)sqrt(( g0 .* k + sig./rho*k.^3 ).*tanh(k*h)); 

%% Finds wave number of suharmonic mode in deep and shallow regions 
options = optimset('Display','off');
Gamma_neutral = @(k,h) sqrt(4./(g0*k*tanh(k*h)).^2.*( (dispEuler(k,h).^2+(2*nu*k.^2).^2 ...
                                    -(omega0/2)^2).^2 + omega0.^2*(2*nu*k.^2).^2));
[kf_deep   ,Gamma_max_deep] = fminsearch(@(k)Gamma_neutral(k,h_deep),1300);
%kf_deep = fsolve(@(k)dispEuler(k,h_deep)-omega0/2,1300,options);
%kf_shallow = fsolve(@(k)dispEuler(k,h_shallow)-dispEuler(kf_deep,h_deep),1300,options); 
[kf_shallow,Gamma_max_shallow] = fminsearch(@(k)Gamma_neutral(k,h_shallow),1300);
lambdaf_deep    =       2*pi/kf_deep;
lambdaf_shallow =       2*pi/kf_shallow;

%% Computes "effective depth" of model
d_deep     = tanh(kf_deep*h_deep) / kf_deep;
d_shallow  = tanh(kf_shallow*h_shallow) / kf_shallow;

%% Drop parameters (dimensional)
num_drops           = 1;                % Number of drops 
drop_type           = 1;                % Type of drop (vector for diff drops)
for n=1:num_drops
    switch drop_type(n)
        case 1 % Fast drop (0.0555 xF/TF)
            drop_radius(n)  = 0.38*10^(-3);  
            theta(n)       = 0.35*2*pi;   %(0.0555 xF/TF) at 0.9 Faraday
            %theta(n)       = 0.322*2*pi; %(0.0575 xF/TF) at 0.95 Faraday
            %theta(n)       = 0.3*2*pi; %(0.0575 xF/TF) at 0.99 Faraday

        case 2 % Medium drop (0.0455 xF/TF)
            drop_radius(n)         = 0.370*10^(-3);
            theta                  = 0.33*2*pi; %Phase, adjusted to match the walker speed

        case 3 % Slow drop 
            drop_radius(n)         = 0.35*10^(-3);
            theta                  = 0.32*2*pi; % (0.0340 xF/TF) at 0.9 Faraday
        case 4 % Very slow drop
            drop_radius(n)         = 0.31*10^(-3);
            theta                  = 0.305*2*pi; % (0.0145 xF/TF) at 0.9 Faraday
    end
end
drop_density        = 949;              % Density of drop (kg/m3)
drop_mass           = 4/3*pi*drop_radius.^3*drop_density; % mass of drop (kg);

%% Air viscosity
mu_air              = 1.8*10^(-5);               % Viscosity of air [kg / (m s)]

%% Choice of scales
TF          = 4*pi/omega0;      % Chosen time scale (Faraday period)
xF          = 2*pi/kf_deep;     % Chosen spatial scale (Farday wavelength)

%% Dimensionless groups
Reynolds    = xF.^2/(TF*nu);        % Reynolds number
nu0         = 1/Reynolds;           % Inverse Reynolds number
Bo          = sig*TF^2/(rho*xF^3);
G           = g0*TF^2/xF;
M           = drop_mass./(rho*xF^3); % 
cf_air      = 6*pi*drop_radius*mu_air*TF/drop_mass;                       % Air drag (vector valued for several drops)
c4          = 0.17; % Coefficient of restitution (Molacek)
cf_impact   = c4*sqrt(rho*drop_radius/sig)*TF*g0;                       % Dissipation during impact

%% Dimensionless depth
h0_deep     = h_deep/xF;
h0_shallow  = h_shallow/xF;

%% Dimensionless FaradayF wavenumber
kf0_deep     = kf_deep*xF;                   
kf0_shallow  = kf_shallow*xF;     

% Dimensionless "effective depth" for wave equation approximation
d0_deep     = d_deep/xF;                    
d0_shallow  = d_shallow/xF;  


%% Gravity
g           = @(t) G*(1 + Gam*cos(4*pi*t-theta));

%% Grid, variable coefficient, and initial data:
Nxy = Nx*Ny;
hx = Lx/Nx; hy = Ly/Ny; Lx = Lx; Ly = Ly;
x = hx*(0:Nx-1)-Lx/2; y = hy*(0:Ny-1)-Ly/2;
[xx,yy] = meshgrid(x,y);
[xxx,yyy] = meshgrid(-Lx/2:0.02:Lx/2-hx,-Ly/2:0.02:Ly/2-hy);


%% Time parameters
dt = dt_desired;
impact_interval = 1;
dt   = impact_interval/ceil(impact_interval/dt);
nsteps_impact = impact_interval/dt;
                      
%% Set the wave-numbers and matrix for multiplication in 2D (Notice i is already included)
kx  =  2*pi*1i/Lx*[0:Nx/2-1 0 -Nx/2+1:-1];
ky  =  2*pi*1i/Ly*[0:Ny/2-1 0 -Ny/2+1:-1];
Kx  = zeros(Ny,Nx); Ky = zeros(Ny,Nx);

for i=1:Ny
    Kx(i,:) = kx;
end
for i=1:Nx
    Ky(:,i) = ky;
end

% Three lines below are to make the program faster
K2 = Kx.^2 + Ky.^2;
abs_K = sqrt(-K2);
dissMatrix = exp(2*nu0*dt*(K2));        % Dissipation operator in FS
dissMatrix_half = exp(2*nu0*dt/2*(K2)); % Dissipation operator in FS
shift1 = mod(-[1:Nx]+1,Nx)+1;
shift2 = mod(-[1:Ny]+1,Ny)+1;
KxiKy  = Kx+1i*Ky;
KxmiKy = Kx-1i*Ky;

phi0 = zeros(size(xx));
eta0 = zeros(size(xx));

%% Plotting options
plotPS = 1;     % Plot power spectrum

%% Store wavefield option
store_wavefield = 1; 
                        
%% Use GPU (CPU typicall faster at low (<512x512) resolution, GPU faster at high resolution)
useGPU =0; 
% Send to GPU
if useGPU ==1
    Kx = gpuArray(Kx);
    Ky = gpuArray(Ky);
    K2 = gpuArray(K2);
    dissMatrix = gpuArray(dissMatrix);
    dissMatrix_half = gpuArray(dissMatrix_half);
    KxiKy = gpuArray(KxiKy);
    KxmiKy = gpuArray(KxmiKy);
    dt = gpuArray(dt);
    nu0 = gpuArray(nu0);
    Bo = gpuArray(Bo);
    shift1 = gpuArray(shift1);
    shift2 = gpuArray(shift2);
end

varList = who;

%Initialiste a structure
p = struct;

%use dynamic fieldnames
for index = 1:numel(varList)
    p.(varList{index}) = eval(varList{index});
end
    
 

 