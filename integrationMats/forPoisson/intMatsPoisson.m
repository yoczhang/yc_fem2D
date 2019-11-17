function [global_intMat, vecRhsF] = ...
    intMatsPoisson(Coeffs_func, func_f, meshInfo, Gaussformulas, basesk)
%
%
%   We let Npoints denote the number of Gauss-Points,
%               Nelems denote the number of the elements of Th,
%               NTbases_trial denote the number of LOCAL trial bases on each K of Th.
%               NTbases_test denote the number of LOCAL test bases on each K of Th.
%
%   input:
%       Coeffs_func, the cell-type, here 
%                   Coeffs_func{1} is the coefficient function,
%                   Coeffs_func{2}, is also the coefficient function, 
%       func_f, the rhs function.
%       meshInfos, the mesh information.
%       formulaGauss2D, the 1d Gauss quadrature formula, size: a matrix, [Npoints x 3],
%               the first column is the x-coordinates of all Gauss-Points,
%               the second column is the y-coordinates of all Gauss-Points,
%               the third is the weights of all Gauss-Points.
%       Tbasesk, the volumn polynomial degree k.
%       Fbasesk, the face polynomial degree k.
%       RTbasesk, the projection polynomial degree k, (== Fbasesk+1).
%
%   output:
%       integrationMat, [(Nelems*NTbases + Nedges*NFbases) x (Nelems*NTbases + Nedges*NFbases)].
%       rhs_fh, [(Nelems*NTbases + Nedges*NFbases) x 1].
%
%
%	YxQian 6/5/2018
%
%   Last modified 7/5/2018
%
%

% mesh setting
Nelems = meshInfo.Nelems;
Nedges = meshInfo.Nedges;
Nnodes = meshInfo.Nnodes;

% coeffs setting
CoeffOne = Coeffs_func{1};
Gaussformula2D = Gaussformulas{1};

% bases setting
%NTbases = nchoosek(basesk+2,basesk);
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

%- global setting
global_intMat = sparse(totalDofs,totalDofs);
vecRhsF = zeros(totalDofs,1);

% the simple build-in function to get the mat
getMatOnElem = @(funcValue, trialValue, testValue, phyGweights) ...
    testValue' * bsxfun(@times, phyGweights.*funcValue, trialValue);
    %> input:
    %>      funcValue, [Npoints x 1], the value the coeffients function at Gauss points.
    %>      trialValue, [Npoints x NTbases_trial], the trial bases values at Gauss points. 
    %>      testValue, [Npoints x NTbases_test], the test bases values at Gauss points. 
    %
    %> output:
    %>       [NTbases x NTbases], the mat of integration ...
    %

for CurrElem = 1:Nelems
    % 1. reference GaussPoints, 
    % 2. different element bases on ref GaussPoints on ii-th element.
    %
    areaElem = meshInfo.areaElem(CurrElem);
    elem = meshInfo.elem{CurrElem};
    vertices = meshInfo.node(elem,:);
    [phyGpoints, ~] = getGaussLocalTri(vertices, Gaussformula2D);
    
    mapMat = [vertices(2,1)-vertices(1,1), vertices(3,1)-vertices(1,1);
        vertices(2,2)-vertices(1,2), vertices(3,2)-vertices(1,2)];
    J_det = 2*areaElem;
    
    [refPb, refPbx, refPby] = femRefBases2D(Gaussformula2D(:,1), Gaussformula2D(:,2), basesk);
        %> refPb, refPbx, refPby, [NTpoints x NTbases].
            
    %--- phy grad basis map to ref grad basis, reference paper: (page:22,
    % section:3.6), 2015 (CMA) FESTUNG: A Matlab/GNU Octave toolbox for the dg.., Part I: Diffusion operator.
    phy2refPbx = 1/J_det*(mapMat(2,2)*refPbx-mapMat(2,1)*refPby); % [NTpoints x NTbases]
    phy2refPby = 1/J_det*(mapMat(1,1)*refPby-mapMat(1,2)*refPbx); % [NTpoints x NTbases]
    
    %%--- the final integration matrix and rhs f.
    %- integration matrix
    valueCoeffOne_2D = CoeffOne(Gaussformula2D(:,1), Gaussformula2D(:,2));
    intMat_Elem = J_det*(...
        getMatOnElem(valueCoeffOne_2D, phy2refPbx, phy2refPbx, Gaussformula2D(:,3)) ...
        + getMatOnElem(valueCoeffOne_2D, phy2refPby, phy2refPby, Gaussformula2D(:,3)) );
        %> [NTbases x NTbases] 
        
%     aa = a(Gaussformula2D(:,1), Gaussformula2D(:,2));
%     intMat_Elem2 = 2*areaElem*(...
%         getMatOnElem(aa.*valueCoeffOne_2D, refPb, refPb, Gaussformula2D(:,3)) );
    
    %- integration rhs f
    value_f = func_f(phyGpoints(:,1), phyGpoints(:,2)); % [NTpoints x 1]
    vec_f = J_det*(refPb'*(Gaussformula2D(:,3).*value_f)); % [NTbases x 1] 
    
    %%--- add local MatElem and rhs_f to global.
    %- add local mat to global
    [Row, Col] = femGetRowColPoisson(CurrElem, meshInfo, basesk);
    
    global_intMat = global_intMat ...
        + sparse(Row(:), Col(:), ...
        intMat_Elem, ...
        totalDofs, totalDofs);
    
    %- add local rhs_f to global
    Row_vec = Row(:,1);
    vecRhsF(Row_vec,1) = vecRhsF(Row_vec,1) + vec_f;
    
end % for CurrElem

end % function matElemCoeffsTrialTest


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

