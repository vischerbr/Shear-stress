function [] = draw_dots(dots, dot_intensities, pars)
% Use this function to make a movie of the dot intensities. Use the trimmed
% dots returned by assign_intensity to cut away gross edge effects.
if(isfield(pars, 'radius'))
    r = pars.radius;
end
if(isfield(pars, 'xmax'))
    r = pars.radius;
end
if(isfield(pars, 'ymax'))
    r = pars.radius;
end
if(isfield(pars, 'sigma'))
    sigma = pars.sigma;
end
if(isfield(pars, 'folder'))
    folder = pars.folder;
end


foldername = ['/home/vischerb/Shear-stress/',folder];
if(~7==exist(foldername))
    mkdir /home/vischerb/Shear-stress/falloff/imgs
end



[xs,ys] = meshgrid(1:1024,1:1024);


for j=1:size(dot_intensities,1)
    j
    img = zeros(1024,1024);
    %sprintf("Now creating image for t=%d", j)
    for l=1:size(dots,1)     
        xc = round(dots(l,1)*1024/pars.xmax);
        yc = round(dots(l,2)*1024/pars.ymax);
        
        img = img + .5*(1+sin(pi*dot_intensities(j,l)))*exp(-((xc-xs).^2+(yc-ys).^2)./(2*sigma^2));
        %draw = rectangle('Position',[px py r*2 r*2],'Curvature',[1,1], 'FaceColor',[1,1,1]*sin(pi/2*dot_intensities(j,i)));
        %daspect([2 2 2])
        
    end
    %frame = getframe;
    img = img./max(max(img));
    img = uint8(round(img*255-1));
    filenametemp = 'test%d.png';
    filename = sprintf(filenametemp, j);
    %[data,map] = frame2im(frame);
    imwrite(img, strcat(foldername,'/imgs/',filename))
end
centroids = dots;
intensity = .5*(1+sin(pi*dot_intensities'));
save([foldername,'intensityrel.mat'], 'intensity')
save([foldername,'location.mat'], 'centroids')
end