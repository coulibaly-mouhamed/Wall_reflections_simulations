function [Fx,Fy] = compute_slope_eta(eta_hat,xi,yi,p)

% Surface gradient
    eta_x_hat = p.Kx.*eta_hat;
    eta_y_hat = p.Ky.*eta_hat;

% Shift fourier spectrum to nearest point
    for n=1:p.num_drops
        ix  = find(p.xx(1,:)>xi(n),1); iy = find(p.yy(:,1)>yi(n),1);
        if isempty(ix); ix=p.Nx; end
        if isempty(iy); iy=p.Ny; end
        shiftx = p.xx(1,ix)-xi(n);
        shifty = p.yy(iy,1)-yi(n);

        eta_x = real(ifft2(exp(-p.Kx.*shiftx-p.Ky.*shifty).*eta_x_hat)); % Shifts solution 
        Fx(n) = eta_x(iy,ix); 
        eta_y = real(ifft2(exp(-p.Kx.*shiftx-p.Ky.*shifty).*eta_y_hat));
        Fy(n) = eta_y(iy,ix);
    end
        
end