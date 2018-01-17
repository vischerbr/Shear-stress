function [frames] = draw_dots(dots, dot_intensities, pars)
% Use this function to make a movie of the dot intensities. Use the trimmed
% dots returned by assign_intensity to cut away gross edge effects.
if(isfield(pars, 'radius'))
    r = pars.radius;
end
foldername = '~/Documents/School/Research/Shear_stress/Test_data/';
if(~7==exist(foldername))
    mkdir Shear_stress/Test_data
end
for j=1:size(dot_intensities,1)
    for i=1:size(dots,1)
        px = dots(i,1) - r;
        py = dots(i,2) - r;
        
        draw = rectangle('Position',[px py r*2 r*2],'Curvature',[1,1], 'FaceColor',[1,1,1]*sin(pi/2*dot_intensities(j,i)));
        daspect([2 2 2])
        hold on;
    end
    frame = getframe;
    
    filenametemp = 'test%d.png';
    filename = sprintf(filenametemp, j);
    %[data,map] = frame2im(frame);
    imwrite(frame.cdata, strcat(foldername,filename))
end

end