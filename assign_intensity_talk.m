function [dot_intensities, dots_trimmed,pacemakers] = assign_intensity_talk(dots, pars)
%ASSIGN_INTENSITY Takes the number of dots and assigns each a sinusoidal
%   intensity profile. No communication between cells. 
%   Detailed explanation goes here

Imax = 1.0; % arbitrary units
dt = .1; % m?
crosstalk_scaling = .05; % arbitrary units
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
I = 2*rand(number_of_dots,1); % arbitrary units

table = make_table(dots,pars);
crosstalktemp = zeros(number_of_dots, number_of_dots);

sprintf("Creating communication matrix.")
for l = 1:number_of_dots
    for q= 1:l
        if q==l
            crosstalktemp(l,q) = self_scaling;%*(-1+2*rand(1));
             
        else    
            crosstalktemp(l,q) = table(l,q)*crosstalk_scaling*rand(1);
        end
    end
end

% If pacemaker cells are enabled, choose a few to make the pacemakers
if(pars.pacemaker)
    pacemakers = zeros(1,15);
    num_pace = 1;
    while num_pace <= size(pacemakers,2)
        testind = randi(number_of_dots);
        if ~ismember(testind,pacemakers)
            pacemakers(num_pace) = testind;
            num_pace = num_pace +1;
        end
    end
    for idx= pacemakers
        crosstalktemp(idx,:) = crosstalktemp(idx,:)*5;
    end
else
    pacemakers = [];
end


% symmetrize the crosstalk
crosstalktemp = crosstalktemp + eye(number_of_dots)*self_scaling;
crosstalk = (crosstalktemp + crosstalktemp');
crosstalk(1:number_of_dots+1:end) = diag(crosstalktemp);
size(crosstalk)


for t=1:tmax
   %sprintf("Calculating intensity for t=%d minutes.", t)
   t
   dI = zeros(number_of_dots,1);
   for p=1:number_of_dots 
       neighbor_table = table - eye(number_of_dots,number_of_dots);
       current_neighbors = neighbor_table(p,:);
       neighbor_intensities = I(find(current_neighbors));
       neighbor_positions = dots(find(current_neighbors),:);
       neighbor_talk = crosstalk(p,find(current_neighbors));
       
       if pars.falloff
            % If we are letting the communication fall off as ~1/r, let's
            % calculate that way
            for f = 1:size(neighbor_intensities)
                neighbor_influence(f) = 1/norm(neighbor_positions(f,:)' - dots(p,:)')*neighbor_talk(f)*sin(pi*(neighbor_intensities(f) - I(p)));
            end
       else
            % else just find the communication normally
            for f = 1:size(neighbor_intensities)
                neighbor_influence(f) = neighbor_talk(f)*sin(pi*(neighbor_intensities(f) - I(p)));
            end
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
sprintf("Trimming border dots away.")
dots_to_trim = zeros(number_of_dots,1);
dots_to_trim = dots(:,1) < R | dots(:,2) < R | dots(:,1) > xmax - R | dots(:,2) > ymax - R;
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

