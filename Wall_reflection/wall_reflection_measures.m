%% set parameters 
L=400;
theta_i = linspace(-90,90,L);
theta_r=  zeros(1,L);
y_a=zeros(1,L);
y_b =zeros(1,L);

%% Run simulation and gather data
n = length(theta_i);

parfor (i= 1:n,40)
    i
    [theta_r(i),y_a(i),y_b(i)] = wall_reflection(theta_i(i)); 
end
save(['reflection_measures','.mat'],...
            'theta_i','theta_r','y_a','y_b');

%% Plot solutions  
%
%plot(theta_i,theta_r*(180/pi));


