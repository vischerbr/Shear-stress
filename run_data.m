pars.nothing = 1;
pars.xmax = 30;
pars.ymax = 30;
pars.radius = .5;
pars.R = 5; %5; % neighbor search radius
pars.ndots = round(.6*pars.xmax*pars.ymax/(pi*pars.radius^2));
pars.tmax = 2000;
pars.dt = .01;
pars.crosstalk = .01; % set higher to account for minimum distance falloff
pars.selftalk = 2.0;
pars.Imax = 1.0;
pars.sigma = 6;
pars.pacemaker = 0;
pars.falloff=1;
pars.folder = 'falloff2';

sprintf("Placing dots.")
dots = place_dots(pars);
pars.ndots = size(dots,1);


%[xs,ys] = meshgrid(1:1024,1:1024);

%sigma = 6;
%intensity = zeros(1024,1024);

%for q =1:ndots
%    xc = round(dots_trimmed(q,1)*1024/pars.xmax);
%    yc = round(dots_trimmed(q,2)*1024/pars.ymax);
    
%    intensity = intensity +exp(-((xc-xs).^2+(yc-ys).^2)./(2*sigma^2));
%end
%intensity_new = intensity./max(max(intensity));
%intensity_new = uint8(round(intensity_new*255-1));
%sprintf("Creating intensity profile.")
[intensities,dots_trimmed,pacemakers] = assign_intensity_talk(dots, pars);
draw_dots(dots_trimmed, intensities, pars)