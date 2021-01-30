M = theta_i<50;
n = sum(M);
L = zeros(1,n);
for i=1:n
    L(i) = (cos(theta_i(i)*pi/180))/(sin(theta_r(i)));
end
plot(theta_i(1:n),L);
xlabel('\theta_i');
ylabel('cos(\theta_i)/sin(\theta_r)');