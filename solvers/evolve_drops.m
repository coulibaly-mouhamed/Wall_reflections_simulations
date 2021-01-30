function [xi, yi, ui, vi] = evolve_drops(xi,yi, ui, vi, p)
% EVOLVES THE DROPS DURING FREE-FLIGHT
% Input:
%   * xi,yi - initial position
%   * ui,vi - initial velocity
%   * p     - structure of needed parameters
% Output:
%   * xi,yi - final position
%   * ui,vi - final velocity

            ui = ui.*exp(-p.cf_air*p.impact_interval);                                          % Speed after free flight
            vi = vi.*exp(-p.cf_air*p.impact_interval);                                          % Speed after free flight
            xi = xi + ui./p.cf_air.*(1-exp(-p.cf_air*p.impact_interval));                       % Position after free flight
            yi = yi + vi./p.cf_air.*(1-exp(-p.cf_air*p.impact_interval));                       % Position after free flight
            xi = mod(xi+p.Lx/2,p.Lx)-p.Lx/2;
            yi = mod(yi+p.Ly/2,p.Ly)-p.Ly/2;
            
end