%% set parameters 
L=301;
theta_i = linspace(0,60,L);
theta_r=  zeros(1,L);
y_a=zeros(1,L);
y_b =zeros(1,L);

%% Run simulation and gather data
n = length(theta_i);

parfor (i= 1:n,40)
    i
    [theta_r_s,y_a_s,y_b_s] = wall_reflection(theta_i(i)); 
    theta_r(i)=theta_r_s;
    y_a(i)=y_a_s;
    y_b(i)= y_b_s;
end
save(['reflection_measures','.mat'],...
            'theta_i','theta_r','y_a','y_b');

%% Plot solutions  
%



