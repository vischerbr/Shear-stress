function [dot_intensities] = assign_intensity_notalk(ndots, tmax)
%ASSIGN_INTENSITY Takes the number of dots and assigns each a sinusoidal
%   intensity profile. No communication between cells.
%   Detailed explanation goes here

Imax = 1.0;
dt = .01;
dot_intensities = zeros(tmax, ndots);
for n= 1:ndots 
    omega = 1+9*rand(1); % rad/sec
    phi = 2*pi*rand(1);
    for t = 1:tmax
        dot_intensities(t,n) = (Imax + Imax*sin(omega*t*dt+phi))/2.0;
        
    end
end

end

