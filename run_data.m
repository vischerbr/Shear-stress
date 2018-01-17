pars.nothing = 1;
pars.xmax = 50;
pars.ymax = 50;
pars.radius = .5;
pars.R = 5; %5; % neighbor search radius
pars.ndots = round(.6*pars.xmax*pars.ymax/(pi*pars.radius^2));
pars.tmax = 100;
pars.dt = .1;
pars.crosstalk = .1;
pars.selftalk = 2.5;
pars.Imax = 1.0;

dots = place_dots(pars);

[intensities,dots_trimmed] = assign_intensity_talk(dots, pars);
draw_dots(dots_trimmed, intensities, pars)