function [Uh, sysInfo] = femPoissonSolve(meshInfo,pde,option)
%
%   We let  Npoints denote the number of Gauss-Points on T,
%               NTbases denote the number of LOCAL bases on each element.
%
%   input:
%       meshInfo, mesh info structure.
%       pde, pde info structure.
%       option, some options about the euqation.
%
%   output:
%       Uh, [Nbases*Nelems x 1], matrix with uh.
%       system, {stiff-matrix,rhs-term}.
%       sysInfo.assembleLocalsolverTime,
%       sysInfo.solveGlobalUhatTime, 
%       sysInfo.solveLocalsolverTime,
%       sysInfo.hdgDiffusionSoverTime.
%
%
%
%	YxQian 6/5/2018
%
%   Last modified 6/5/2018
%


if ~exist('option','var')
    error('there is no option in femPoissonSolve.m')
end

basesk = basesType2degreek(option.basestype);
Nelems = meshInfo.Nelems;
Nedges = meshInfo.Nedges;
Nnodes = meshInfo.Nnodes;
V_Ndofs = Nnodes;
if basesk >= 2
    oneFdofs = nchoosek((basesk-2)+1,(basesk-2));
    F_Ndofs = Nedges * oneFdofs;
else 
    F_Ndofs = 0;
end
if basesk >= 3
    oneTdofs = nchoosek((basesk-3)+2,(basesk-3));
    T_Ndofs = Nelems * oneTdofs;
else
    T_Ndofs = 0;
end
totalDofs = V_Ndofs + F_Ndofs + T_Ndofs;


%%  generate the Gaussformulas
% 2d 
[Q1,Q2,W] = quadRule2D(2*basesk+1); % Gauss-Points on Ref elem
Gaussformula2D = [Q1', Q2', W']; 
    %> In order to facilitate the calculation, we structure the 'formula' matrix,
    %> [Npoints x 3],
    %> the first column is the x-coordinates of all Gauss-Points,
    %> the second column is the y-coordinates of all Gauss-Points,
    %> the third is the weights of all Gauss-Points.

% 1d
% [q,w] =quadRule1D(2*trial_k+1); % Gauss-Points on [0,1]
[q,w] =quadRule1D(12); % Gauss-Points on [0,1]
q = 2*q-1; w = 2*w; % Gauss-Points [0,1] --> [-1,1].
Gaussformula1D = [q', w'];
    %> [Npoints x 2],
    %> the first column is the 1D coordinates of all Gauss-Points,
    %> the second is the weights of all Gauss-Points.

% get the uniform Gaussformulas, is a CELL structure data. 
Gaussformulas = {Gaussformula2D,... % 2d quad formula.
    Gaussformula1D}; % 1d quad formula.


%% edge information
meshInfo = getPoissonBoundaryInfo(meshInfo);

%% assemeble matrix and rhs
tic; %<<<<<<<<<<<<<< tic2
Coeffs_func=cell(3,1);
Coeffs_func{1} = pde.funcOne;
Coeffs_func{2} = pde.funcOne;
Coeffs_func{3} = pde.funcOne;

%%------- assemble the matrix ------%%
[sysM, rhs_fh] = ...
    intMatsPoisson(Coeffs_func, pde.f, meshInfo, Gaussformulas, basesk);

%%------- treat the boundary condition ------%%
DirichletDofs = femFindBoundaryDofs(meshInfo, basesk);

if ~isempty(DirichletDofs)
    uD = femTreatDirichletBCPoisson(pde.gD, meshInfo, basesk);
end
if ~isempty(meshInfo.NeumannEdgeIndex)
    uN = femTreatNeumannBCPoisson(pde, meshInfo, Gaussformula1D, basesk);
    rhs_fh = rhs_fh + uN;
end

%--- assignment the sysInfo
sysInfo.AssembleTime = toc; %<<<<<<<<<<<<<< corresponding tic2
disp(['assemble time: ',num2str(sysInfo.AssembleTime)])


%% solve the system
tic; %<<<<<<<<<<<<<< tic3
%%------- first, get the solution of dofs on edges ------%%
Uh = zeros(totalDofs, 1);
Uh(DirichletDofs) = uD;

freeDofs = (1:totalDofs)';
freeDofs(DirichletDofs) = [];

rhs_fh = rhs_fh - sysM(:,DirichletDofs) * Uh(DirichletDofs);
Uh(freeDofs) = sysM(freeDofs,freeDofs)\rhs_fh(freeDofs);

sysInfo.SolveSystemTime = toc; %<<<<<<<<<<<<<< corresponding tic3
disp(['solve the system matrix time: ',num2str(sysInfo.SolveSystemTime)])

%% assignment the sysInfo
sysInfo.Ndofs = totalDofs;
sysInfo.basesk = basesk;
sysInfo.Gaussformulas = Gaussformulas; 

end % function


%% ------- sub functions -------
%
%
%
%%------- getStokesBoundaryInfo, sub function 1 -------%%
function meshInfo = getPoissonBoundaryInfo(meshInfo)
%
%
%
%
%
%	YxQian 6/5/2018
%
%   Last modified 6/5/2018
%


%--- the domain is [0,1]x[0,1]
bdEdgeIndex = meshInfo.bdEdgeIndex; % here, all the bdEdge 

%--- Dirichlet boundary
DirichletEdgeIndex_1 = bdEdgeIndex( abs(meshInfo.baryEdge(bdEdgeIndex,1)-0) <= 5e-8 ); % left bd
DirichletEdgeIndex_2 = bdEdgeIndex( abs(meshInfo.baryEdge(bdEdgeIndex,1)-1) <= 5e-8 ); % right bd
DirichletEdgeIndex_3 = bdEdgeIndex( abs(meshInfo.baryEdge(bdEdgeIndex,2)-1) <= 5e-8 ); % top bd
meshInfo.DirichletEdgeIndex = ...
    union([DirichletEdgeIndex_1;DirichletEdgeIndex_2],DirichletEdgeIndex_3);
%meshInfo.DirichletEdgeIndex = bdEdgeIndex;

%--- Neumann boundary
NeumannEdgeIndex_1 = bdEdgeIndex( abs(meshInfo.baryEdge(bdEdgeIndex,2)-0) <= 5e-8 ); % bottom bd
%meshInfo.NeumannEdgeIndex = union(NeumannEdgeIndex_1,NeumannEdgeIndex_2);
meshInfo.NeumannEdgeIndex = NeumannEdgeIndex_1;
%meshInfo.NeumannEdgeIndex = [];
end % function 


function DirichletDofs = femFindBoundaryDofs(meshInfo, basesk)
%
%
%
%

%--- Dirichlet dofs
DirEdgeIndex = meshInfo.DirichletEdgeIndex;
DirVdofs = meshInfo.edge(DirEdgeIndex,:);
DirVdofs = unique(DirVdofs(:));

if basesk >= 2
    oneFdofs = nchoosek((basesk-2)+1,(basesk-2));
    DirFdofs = meshInfo.Nnodes + (DirEdgeIndex-1) * oneFdofs; % (First d.o.f.) of each Dirichlet face.
    DirFdofs = bsxfun(@plus,DirFdofs,1:oneFdofs); % [Ndir x oneFdofs] (d.o.f. for each Dirichlet face) 
    DirFdofs = DirFdofs';
    DirFdofs = DirFdofs(:);
    % %DirFdofs = reshape(DirFdofs', oneFdofs*length(DirEdgeIndex), 1);
else 
    DirFdofs = [];
end

DirichletDofs = [DirVdofs;
    DirFdofs];

end % function