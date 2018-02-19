function [ neighbor_table ] = make_table(dots,pars)
%make_table Dirtily constructs neighbor table for created dots. Radius is
%currently hardcoded, I'm not a software engineer.
%   Detailed explanation goes here

% Construct output matrix
number_of_dots = size(dots,1);
neighbor_table = zeros(number_of_dots,number_of_dots);


R = pars.R;

% Use range search to create list of nearest neighbors
if pars.falloff
    neighbor_search = rangesearch(dots,dots,3*R);
else
    neighbor_search = rangesearch(dots,dots,R);
end
% Loop over list of neighbors and construct logical neighbor table

for m = 1:length(neighbor_search)
    m_neighbors = cell2mat(neighbor_search(m));
    for n =1:number_of_dots
        if ismember(n, m_neighbors)
            neighbor_table(m,n)=1;
        end 
    end
end
%neighbor_table = neighbor_table - eye(number_of_dots,number_of_dots); % remove the diagonal
%diagonal currently remaining to let crosstalk matrix handle self talk
end

