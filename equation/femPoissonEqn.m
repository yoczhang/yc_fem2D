function [sysErr,sysTime] = femPoissonEqn(pde,option)
%
%   femPoissonEqn solve the poisson equation by FEM methods
% 
%
%	YxQian 6/5/2018
%
%   Last modified 6/5/2018
%

%% Default setting of mesh and pde data
if ~exist('option','var')
    error('there is no option in femPoissonEqn.m')
end
if ~exist('pde','var')
    pde = poissonData(0,5); % default data
end


%% Parameters
maxIt = option.maxIt;


%% Initialize err
sysErr = zeros(maxIt,1); % not used
sysTime = zeros(maxIt,1); % not used

uh_L2_error = zeros(maxIt,1); uh_H1_error = zeros(maxIt,1);
uhL2rate = zeros(maxIt,1); uhH1rate = zeros(maxIt,1);
h = zeros(maxIt,1);
% Nnodes = zeros(maxIt,1);
% Ndof = zeros(maxIt,1);


%% set path and file name
% load('setpath_pwd.mat')
% femfunc_name = 'femPoisson_';

%>>>>>>>>>  creat log file >>>>>>>>
% date = datestr(now,31); 
%     %> capture the now time, 31 stands for the scheme: 2000-03-01 15:45:17
% logFilename=[setpath_pwd,'/logs/',femfunc_name,date,'_log.txt'];
% diary(logFilename);
% diary on; % begin diary
%<<<<<<<<<<<<<<<<<<<<<<<<<<

disp('********************  fem Poisson  ********************')
disp('------------------------------------------------------')
disp('poisson Equation:')
disp(['   K_11 = ',num2str( pde.K_11)])
disp(['   K_22 = ',num2str( pde.K_22)])
disp(['   u = ',func2str(pde.u)])
disp('FEM parameters: ')
disp(['   bases: ', option.basestype])
disp(['   maxIt: ', num2str(maxIt)])

%% FEM methods
disp('------------------------------------------------------')
disp('Solving hhoPoisson: ')
disp('------------------------------------------------------')
format long e
for n = 1:maxIt
    disp(' ')
    n_str = num2str(n);
    n_th_cycle = ['the ',n_str,'-th cycle(i.e. the ', n_str, '-th mesh)'];
    disp(n_th_cycle)
    
    %% get mesh information
    %---- mesh 
	%-------------------- Many mesh ---------------------
	meshname = 'tri';
	Dmeshname = ['Dmesh_',meshname,'_[0,1]x[0,1]_',num2str(2^(n+1))];
	load(Dmeshname); 
	%-----------------------------------------------------
    
    %--- mesh 2
    %-------------------- Tri mesh ---------------------
%     g_mesh = domainSquare(1/2^(n+1));
%     node = g_mesh.coordV;
%     elem = g_mesh.V0T;
    %-----------------------------------------------------

    meshInfo = polyMeshAuxStructure(node, elem);
    %patchPlotMesh(node, elem);
    %plotPolyMsh(meshInfo)

    %% solve equations
    solve_t0 = cputime;
    [Uh, sysInfo]= femPoissonSolve(meshInfo,pde,option);
    sysInfo.SoverTime = cputime - solve_t0;
    disp(['the solve-func time: ',num2str(sysInfo.SoverTime)])
    
    basesk = sysInfo.basesk;
    Gaussformula2D = sysInfo.Gaussformulas{1};
    %% compute the err
    err_t0 = cputime;
    [uh_L2_error(n), uh_H1_error(n)] = femL2H1error(pde.u,pde.ux,pde.uy,Uh,meshInfo,Gaussformula2D,basesk);
    sysInfo.ErrTime = cputime - err_t0;
    disp(['compute Err time: ',num2str(sysInfo.ErrTime)])
    disp(['uh_L2_error = ',num2str(uh_L2_error(n))])
    disp(['uh_H1_error = ',num2str(uh_H1_error(n))])
    
    %% other options
    h(n) = 1./(sqrt(meshInfo.Nnodes)-1);

    disp('------------------------------------------------------')
    disp('--- show solution ----')
    showsolution(meshInfo.node, cell2mat(meshInfo.elem), Uh);
end % for n

disp('------------------------------------------------------')
%%
uhL2rate(2:maxIt) = log(uh_L2_error(1:maxIt-1)./uh_L2_error(2:maxIt)) ...
    ./log(h(1:maxIt-1)./h(2:maxIt));
uhH1rate(2:maxIt) =  log(uh_H1_error(1:maxIt-1)./uh_H1_error(2:maxIt)) ...
    ./log(h(1:maxIt-1)./h(2:maxIt));

disp('Table: Error')
colname = {'h  ', '   ||u-U_h||_0 ','   ||u-U_h||_1'};
disptable(colname, h,'%0.2e', uh_L2_error,'%0.5e', ...
    uh_H1_error,'%0.5e');

disp('Table: Error rate')
colname = {'h', '   UhL2rate', '   UhH1rate'};
disptable(colname, h,'%0.2e', uhL2rate,'%0.4f', ...
    uhH1rate,'%0.4f');
disp('------------------------------------------------------')

%>>>>>>>>>  close the diary >>>>>>
% diary off; % close the diary
%<<<<<<<<<<<<<<<<<<<<<<<<<<
end % function dgPoissonEqn
