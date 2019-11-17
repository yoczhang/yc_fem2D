% Generates a mesh for the unit square.
%

%
function g = domainSquare(h)
dim = ceil(1/h); % number of edges per side of the unit square
h = 1/dim;
%% Build coordV.
[X, Y] = meshgrid(0:h:1);
Xlist = reshape(X, length(X)^2, 1);  Ylist = reshape(Y, length(X)^2, 1);
coordV = [Xlist, Ylist];
%% Build V0T.
pat1 = [1,dim+2,2]; % pattern of "lower-left" triangles
V0T1 = repmat(pat1, dim*(dim+1), 1) + repmat((0:dim*(dim+1)-1)', 1, 3);
V0T1(dim+1 : dim+1 : dim*(dim+1), :) = [];
pat2 = [dim+2,dim+3,2];
V0T2 = repmat(pat2, dim*(dim+1), 1) + repmat((0:dim*(dim+1)-1)', 1, 3);
V0T2(dim+1 : dim+1 : dim*(dim+1), :) = [];
%% Generate grid data and boundary IDs
g = generateGridData(coordV, [V0T1; V0T2]);
g.idE = zeros(g.numE, 1);
g.idE(g.baryE(:, 2) == 0) = 1; % south
g.idE(g.baryE(:, 1) == 1) = 1; % east
g.idE(g.baryE(:, 2) == 1) = 1; % north
g.idE(g.baryE(:, 1) == 0) = 1; % west
g.idE0T = g.idE(g.E0T); % local edge IDs
end % function
