
theta_r =theta_r*180/pi;
theta_r = theta_r -fix(theta_r/90)*90;
plot(theta_i,theta_r);
xlabel('\theta_i');
ylabel('\theta_r');