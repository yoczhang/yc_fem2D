function [L2err, H1err] = femL2H1error(U,Ux,Uy,Uh,meshInfo,Gaussformula2D,basesk)
%
%   We let Npoints denote the number of Gauss-Points,
%               Nelems denote the number of the elements of Th,
%               NTbases denote the number of LOCAL bases on each K of Th.
%
%   input:
%       U, vectorized function of two variables (x,y), i.e. the true solutions.
%       Ux, vectorized function of two variables (x,y), i.e. the true solutions of partial derivate x.
%       Uy, vectorized function of two variables (x,y), i.e. the true solutions of partial derivate y.
%       Uh, discontinuous Pk function, [totalDofs x 1]
%       meshInfo, the mesh information.
%       formulaGauss2D, the 2d Gauss quadrature formula, size: a matrix, Npoints x 3,
%               the first column is the x-coordinates of all Gauss-Points,
%               the second column is the y-coordinates of all Gauss-Points,
%               the third is the weights of all Gauss-Points.
%       basesk, polynomial degree.
%
%   output:
%       L2err, \| P-Ph \|_{L^2(Th)}, a scalar.
%       H2err, \| P-Ph \|_{H^1(Th)}, a scalar.
%
%
%	YxQian 7/5/2018
%
%   Last modified 7/5/2018
%

% mesh information: interior edges 
Nelems = meshInfo.Nelems;

L2err = 0;
H1err = 0;

for CurrElem = 1:Nelems
    %% Part I, get the information about 
    % 1. physical GaussPoints, 
    % 2. different element bases on phy GaussPoints on ii-th element.
    %
    %>>-- Begin Part I -------------------------------- DONOT MODIFY ------------------------------
    areaElem = meshInfo.areaElem(CurrElem);
    elem = meshInfo.elem{CurrElem};
    vertices = meshInfo.node(elem,:);
    [phyGpoints, phyGweights] = getGaussLocalTri(vertices, Gaussformula2D);
    
    mapMat = [vertices(2,1)-vertices(1,1), vertices(3,1)-vertices(1,1);
        vertices(2,2)-vertices(1,2), vertices(3,2)-vertices(1,2)];
    J_det = 2*areaElem;
    
    [refPb, refPbx, refPby] = femRefBases2D(Gaussformula2D(:,1), Gaussformula2D(:,2), basesk);
        %> refPb, refPbx, refPby, [NTpoints x NTbases].
            
    %--- phy grad basis map to ref grad basis, reference paper: (page:22,
    % section:3.6), 2015 (CMA) FESTUNG: A Matlab/GNU Octave toolbox for the dg.., Part I: Diffusion operator.
    phy2refPbx = 1/J_det*(mapMat(2,2)*refPbx-mapMat(2,1)*refPby); % [NTpoints x NTbases]
    phy2refPby = 1/J_det*(mapMat(1,1)*refPby-mapMat(1,2)*refPbx); % [NTpoints x NTbases]
    
    %--- get local dofs
    [Row, ~] = femGetRowColPoisson(CurrElem, meshInfo, basesk);
    localDofs = Uh(Row(:,1)); % [NTbases x 1]
    
    %--- value of true solution at Gauss-points
    value_U = U(phyGpoints(:,1), phyGpoints(:,2)); % [NTpoints x 1]
    value_Ux = Ux(phyGpoints(:,1), phyGpoints(:,2)); % [NTpoints x 1]
    value_Uy = Uy(phyGpoints(:,1), phyGpoints(:,2)); % [NTpoints x 1]
    
    %--- L2err, H1err
    L2err = L2err + phyGweights'*(value_U - refPb*localDofs).^2;
    H1err = H1err + phyGweights'*(value_Ux - phy2refPbx*localDofs).^2 ...
        +  phyGweights'*(value_Uy - phy2refPby*localDofs).^2;

end % for CurrElem
H1err = sqrt(H1err);
L2err = sqrt(L2err);

end % function dgL2H1Error




%% ------------------------------------------ sub function -------------------------------------------------- %%
%--------------------------------------------------------------------------------------------------------------------%
%%>> -- Begin sub function 1 -------------------------------------------------------------
function [phyGpoints, phyGweights] = getGaussLocalTri(coordTri_nt, formulaGauss2D)
%
%   output:
%       phyGpoints, [Npoints x 2]
%       phyGweights, [Npoints x 1]
x1=coordTri_nt(1,1);
y1=coordTri_nt(1,2);
x2=coordTri_nt(2,1);
y2=coordTri_nt(2,2);
x3=coordTri_nt(3,1);
y3=coordTri_nt(3,2);
JacobiTri=abs((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1));

phyGweights = JacobiTri * formulaGauss2D(:,3);
phyGpoints(:,1)=x1+(x2-x1)*formulaGauss2D(:,1)+(x3-x1)*formulaGauss2D(:,2);
phyGpoints(:,2)=y1+(y2-y1)*formulaGauss2D(:,1)+(y3-y1)*formulaGauss2D(:,2);
end % function getGaussLocalTri
%%<<-- End sub function 2 ---------------------------------------------------------------
