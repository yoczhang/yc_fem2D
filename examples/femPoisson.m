%
%   RATE OF CONVERGENCE OF FEM METHOD FOR THE DIFFUSION EQNs IN 2D
%
%   This example is to show the rate of convergence of FEM
%	approximation of the Diffusion equation on the unit square:
%
%       -div(grad u) = f in \Omega,
%                        u = gD  on \Gamma_D,
%     grad u \cdot n = (gN1,gN2)\cdot n on \Gamma_N.
%
%
%
%
%	YxQian 6/5/2018
%
%   Last modified 6/5/2018
%


clc
clearvars; close all;
%% Setting
% program setting
option.L0 = 1; % to uniformrefine the mesh and generate an initial mesh 
option.maxIt = 4; % to define the maximum circulation, to generate the Rate.
option.printlevel = 1;
option.plotflag = 1;

% bases element setting
option.basestype = 'P2'; % may take 'k-1','k','k+1';

%% pde setting
pde = poissonData(0,12);
    % input: (coeff_case, u_case)

%% equation
femPoissonEqn(pde,option);

%%-- end function