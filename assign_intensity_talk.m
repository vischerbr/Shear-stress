function [dot_intensities, dots_trimmed] = assign_intensity_talk(dots, pars)
%ASSIGN_INTENSITY Takes the number of dots and assigns each a sinusoidal
%   intensity profile. No communication between cells. 
%   Detailed explanation goes here

Imax = 1.0; % arbitrary units
dt = .1; % s
crosstalk_scaling = .1; % arbitrary units
self_scaling = 2.5;

if(isfield(pars, 'xmax'))
    xmax = pars.xmax;
end
if(isfield(pars, 'ymax'))
    ymax = pars.ymax;
end
if(isfield(pars, 'dt'))
    dt = pars.dt;
end
if(isfield(pars, 'Imax'))
    Imax = pars.Imax;
end
if(isfield(pars, 'crosstalk'))
    crosstalk_scaling = pars.crosstalk;
end
if(isfield(pars, 'selftalk'))
    self_scaling = pars.selftalk;
end
if(isfield(pars, 'tmax'))
    tmax = pars.tmax;
end
if(isfield(pars, 'R'))
    R = pars.R;
end


number_of_dots = size(dots,1);



dot_intensities = zeros(tmax, number_of_dots);
I = rand(number_of_dots,1); % arbitrary units

table = make_table(dots);
crosstalktemp = zeros(number_of_dots, number_of_dots);

for l = 1:number_of_dots
    for q= 1:l
        if q==l
            crosstalktemp(l,q) = self_scaling*(-1+2*rand(1));
        else    
            crosstalktemp(l,q) = table(l,q)*crosstalk_scaling*rand(1);
        end
    end
end

% symmetrize the crosstalk

crosstalk = (crosstalktemp + crosstalktemp');
crosstalk(1:number_of_dots+1:end) = diag(crosstalktemp);
size(crosstalk)
for t=1:tmax
   dI = zeros(number_of_dots,1);
   for p=1:number_of_dots 
       current_neighbors = table(p,:);
       neighbor_intensities = I(find(current_neighbors));
       neighbor_talk = crosstalk(p,find(current_neighbors));
       
       for f = 1:size(neighbor_intensities)
           neighbor_influence(f) = neighbor_talk(f)*sin(pi*(neighbor_intensities(f) - I(p)));
       end
           
       %arrayfun(@(x) sin(pi*(x- I(p))), neighbor_intensities)
       dI(p) = dt*(crosstalk(p,p) + sum(neighbor_influence));
   end    
   dot_intensities(t,:) = I;
   I = mod(I + dI,2.0);
end


% Pare down the cells to diminish edge effects - some reside, but this is a
% good bandage fix. Note this enforces that xmax, ymax > R, but you should
% be doing that anyway.
dots_to_trim = zeros(number_of_dots,1);
dots_to_trim = dots(:,1) < R | dots(:,2) < R | dots(:,1) > xmax - R | dots(:,2) > ymax - R
dots(dots_to_trim,:) = [];
dots_trimmed = dots;
dot_intensities(:,dots_to_trim) = [];
% 
% for n= 1:number_of_dots 
%     omega = 1+9*rand(1); % rad/sec
%     phi = 2*pi*rand(1);
%     for t = 1:tmax
%         dot_intensities(t,n) = (Imax + Imax*sin(omega*t*dt+phi))/2.0;
%         
%     end
% end
% 
% end
end

