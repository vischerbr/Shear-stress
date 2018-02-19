function [dots] = place_dots(pars)
%PLACE_DOTS Places dots randomly in an area, used to generate arbitrary
% cell intensity profiles. Stolen librally from stack exchange.
%   Please write me, I explain how this is used! 

if(isfield(pars, 'xmax'))
    xmax = pars.xmax;
end

if(isfield(pars, 'ymax'))
    ymax = pars.ymax;
end

if(isfield(pars, 'ndots'))
    ndots = pars.ndots;
end

if(isfield(pars, 'radius'))
    radius = pars.radius;
end

dots = zeros(ndots ,2);

for i=1:ndots
    %Flag which holds true whenever a new dot was found
    newdotFound = false;
    j=0;
    %loop iteration which runs until finding a dot which doesnt intersect with previous ones
    while ~newdotFound
        if j>10000
            break;
        end
        x = 0 + (xmax)*rand(1);
        y = 0 + (ymax)*rand(1);
        %calculates distances from previous drawn dots
        prevdotsX = dots(1:i-1,1);
        prevdotsY = dots(1:i-1,2);
        distFromPrevdots = ((prevdotsX-x).^2+(prevdotsY-y).^2).^0.5;

        %if the distance is not to small - adds the new dot to the list
        if i==1 || sum(distFromPrevdots<=2*radius)==0
            newdotFound = true;
            dots(i,:) = [x y];
        end
        j=j+1;
    end
    if j>10000
            break;
    end
end

dots = dots(any(dots,2),:);
sprintf("Placed %d dots.", size(dots, 1))
end

