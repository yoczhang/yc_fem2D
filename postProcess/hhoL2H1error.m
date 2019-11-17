function [L2err, H1err] = hhoL2H1error(P,Ph_T,Ph_F,meshInfo,Gaussformulas,Tbasesk,Fbasesk)
%
%   %------------------------------------------------------------------
%       [5/1/2018]
%       In this file, we compute the L2 err and H1 err,
%       which, the L2 err is just on the interior of element T,
%       but, the H1 err is combine the Face and Element err.
%   %-------------------------------------------------------------------
%
%   We let Npoints denote the number of Gauss-Points,
%               Nelems denote the number of the elements of Th,
%               NTbases denote the number of LOCAL bases on each K of Th.
%
%   input:
%       P, vectorized function of two variables (x,y), i.e. the true solutions.
%       Ph, discontinuous Pk function, [NTbases*Nelems x 1]
%       meshInfo, the mesh information.
%       formulaGauss2D, the 2d Gauss quadrature formula, size: a matrix, Npoints x 3,
%               the first column is the x-coordinates of all Gauss-Points,
%               the second column is the y-coordinates of all Gauss-Points,
%               the third is the weights of all Gauss-Points.
%       basesType_trial, polynomial degree.
%
%   output:
%       L2err, \| P-Ph \|_{L^2(Th)}, a scalar.
%       H2err, \| P-Ph \|_{H^1(Th)}, a scalar.
%
%
%   YcZhang 5/1/2018
%
%   Last modified 5/1/2018
%

Gaussformula2D = Gaussformulas{1};
Gaussformula1D = Gaussformulas{2};

% mesh information: interior edges 
Nelems = meshInfo.Nelems;
NTbases = (Tbasesk+1)*(Tbasesk+2)/2;
NFbases = nchoosek(Tbasesk+1,Tbasesk);

L2err = 0;
H1err = 0;

for CurrElem = 1:Nelems
    %% Part I, get the information about 
    % 1. physical GaussPoints, 
    % 2. different element bases on phy GaussPoints on ii-th element.
    %
    %>>-- Begin Part I -------------------------------- DONOT MODIFY ------------------------------
    singleElem = meshInfo.elem{CurrElem,1}; % get the ii-th elem from the cell-type structure: elem.
    singleNE = size(singleElem,2); % get the number of edges on singleElem.
    
    concavePointElem_ii = meshInfo.concavePointElem(CurrElem); % [1 x 1], the concavePointElem of ii-the elem.
    coordv = meshInfo.node(singleElem,:); % [singleNE x 2], the vertices (x-coord, y-coord) of element.
    centroidElem = meshInfo.centroidElem(CurrElem,:); % [1 x 2], the centroid(xing xin) (x-coord, y-coord) of element.
    
    coordTri0Elem = getcoordTri0Elem(singleNE, concavePointElem_ii, centroidElem, coordv); % if singleNE==3, coordTri0Elem, [3 x 2]. else coordTri0Elem, [3*singleNE x 2].
    
    elem_xT = centroidElem(1);  elem_yT = centroidElem(2); 
    elem_hT = meshInfo.hElem(CurrElem); 
    
    PhOnElem = Ph_T((CurrElem-1)*NTbases+1 : CurrElem*NTbases,1); % [NTbases_trial x 1]
    %<<-- End Part I ---------------------------------------------------------------------------------------
    
    mat_temp = zeros(NTbases,NTbases);
    valueP_temp = zeros(NTbases,1);
    for nt = 1:(size(coordTri0Elem,1)/3)
        coordTri_nt = coordTri0Elem((nt-1)*3+1:nt*3,:); % [3 x 2], (x-coorc,y-coord) of coordTri_nt.
        [phyGpoints, phyGweights] = getGaussLocalTri(coordTri_nt, Gaussformula2D);

        %--- get the bases values on quad points
        [uT_Pb, ~, ~] = localBases2D(elem_xT, elem_yT, elem_hT, phyGpoints(:,1), phyGpoints(:,2), Tbasesk);
        wuTPb = bsxfun(@times, phyGweights, uT_Pb); % trialPb, trialPbx, trialPby, [Npoints x NTbases_trial]      
            
        %--- the funcValue may be chosen by case.
        valueP_T = P(phyGpoints(:,1), phyGpoints(:,2)); % [Npoints x 1]
        
        mat_temp = mat_temp + (wuTPb'*uT_Pb);
        valueP_temp = valueP_temp + wuTPb'*valueP_T;
    end % for nt
    
    %--- get the of exact solution
    Pt_Th = mat_temp\valueP_temp; % the L2 projection of exact solution on Elem.
    
    for nt = 1:(size(coordTri0Elem,1)/3)
        coordTri_nt = coordTri0Elem((nt-1)*3+1:nt*3,:); % [3 x 2], (x-coorc,y-coord) of coordTri_nt.
        [phyGpoints, phyGweights] = getGaussLocalTri(coordTri_nt, Gaussformula2D);

        % get the bases values on quad points
        [uT_Pb, uT_Pbx, uT_Pby] = localBases2D(elem_xT, elem_yT, elem_hT, ...
            phyGpoints(:,1), phyGpoints(:,2), Tbasesk);% trialPb, trialPbx, trialPby, [Npoints x NTbases_trial]      
        
        %--- L2err
        L2err = L2err + phyGweights'*(uT_Pb*Pt_Th - uT_Pb*PhOnElem).^2;
        
        %--- H1err, part I.
        H1err = H1err + phyGweights'*(uT_Pbx*Pt_Th-uT_Pbx*PhOnElem).^2 ...
            + phyGweights'*(uT_Pby*Pt_Th-uT_Pby*PhOnElem).^2;

    end % for nt
    
    for CurrEdge = 1:singleNE
        %--- edge setting
        edgeIndex = meshInfo.elem2edge{CurrElem}(CurrEdge);
        ePoint1 = meshInfo.node(meshInfo.edge(edgeIndex,1),:); % [1 x 2] the (x-coordinate,y-coordinate) of ePoint1
        ePoint2 = meshInfo.node(meshInfo.edge(edgeIndex,2),:); % [1 x 2] the (x-coordinate,y-coordinate) of ePoint2
        edge_xE = (ePoint1(1) + ePoint2(1))/2;
        %edge_yE = (ePoint1(2) + ePoint2(2))/2;
        edge_hE = meshInfo.areaEdge(edgeIndex);
        
        PhOnFace = Ph_F((NFbases*(edgeIndex-1)+1):NFbases*edgeIndex);
        
        %--- 1D Gauss Points
        phyGpoints1DX = (ePoint1(1)-ePoint2(1))/2*Gaussformula1D(:,1) + (ePoint1(1)+ePoint2(1))/2; 
        phyGpoints1DY = (ePoint1(2)-ePoint2(2))/2*Gaussformula1D(:,1) + (ePoint1(2)+ePoint2(2))/2; 
        phyGweights1D = edge_hE*Gaussformula1D(:,2)/2;
        
        %--- some values at 1D GaussPoints
        [uT_Pb, ~, ~] = localBases2D(elem_xT, elem_yT, elem_hT, phyGpoints1DX, phyGpoints1DY, Tbasesk);
            %> T_Pb, T_Pbx, T_Pby, [NFpoints x uNTbases].
        [uF_Pb, ~] = localBases1D(edge_xE, edge_hE, Gaussformula1D(:,1), Fbasesk);
            %> F_Pb, F_Pbx, F_Pby, [NFpoints x uNFbases].
            
        wuFPb = bsxfun(@times, phyGweights1D, uF_Pb); % [NFpoints x NFbases]
        valueP_F = P(phyGpoints1DX, phyGpoints1DY); % [NFpoints x 1]
        Pt_Fh = (wuFPb'*uF_Pb)\(wuFPb'*valueP_F); % the exact solution P, projection on the face.
        
        %--- H1err, part II.
        H1err = H1err + ...
            edge_hE^(-1)*phyGweights1D'*((uT_Pb*PhOnElem-uT_Pb*Pt_Th) ...
            - (uF_Pb*PhOnFace-uF_Pb*Pt_Fh)).^2;
    end

end % for CurrElem
% H1err = sqrt(L2err+H1err)/sqrt(projection_Pthx+projection_Pthy);
% L2err = sqrt(L2err)/sqrt(projection_Pth);
H1err = sqrt(H1err);
L2err = sqrt(L2err);

end % function dgL2H1Error




%% ------------------------------------------ sub function -------------------------------------------------- %%
%--------------------------------------------------------------------------------------------------------------------%
%%>>-- Begin sub function 1 ---------------------------------------------------------
function coordTri0Elem = getcoordTri0Elem(singleNE, beginP_n, baryElem, coordv)
%
%
%   input:
%       beginP_n, the n-th point of singleElem, and from the n-th point to
%       construct the little triangles.
%

m = @(x) mod(x,singleNE)+(x==singleNE)*singleNE;
n = beginP_n;

if n == 0
    if singleNE == 3
        coordTri0Elem = coordv; % [3 x 2]
    elseif singleNE == 4
        coordTri0Elem = ...
            [baryElem; coordv(1,:); coordv(2,:);
            baryElem; coordv(2,:); coordv(3,:);
            baryElem; coordv(3,:); coordv(4,:);
            baryElem; coordv(4,:); coordv(1,:)]; % [3*singleNE x 2]
    elseif singleNE == 5
        coordTri0Elem = ...
            [baryElem; coordv(1,:); coordv(2,:);
            baryElem; coordv(2,:); coordv(3,:);
            baryElem; coordv(3,:); coordv(4,:);
            baryElem; coordv(4,:); coordv(5,:);
            baryElem; coordv(5,:); coordv(1,:)];
    elseif singleNE == 6
        coordTri0Elem = ...
            [baryElem; coordv(1,:); coordv(2,:);
            baryElem; coordv(2,:); coordv(3,:);
            baryElem; coordv(3,:); coordv(4,:);
            baryElem; coordv(4,:); coordv(5,:);
            baryElem; coordv(5,:); coordv(6,:);
            baryElem; coordv(6,:); coordv(1,:)]; % [3*singleNE x 2]
    elseif singleNE == 7
        coordTri0Elem = ...
            [baryElem; coordv(1,:); coordv(2,:);
            baryElem; coordv(2,:); coordv(3,:);
            baryElem; coordv(3,:); coordv(4,:);
            baryElem; coordv(4,:); coordv(5,:);
            baryElem; coordv(5,:); coordv(6,:);
            baryElem; coordv(6,:); coordv(7,:);
            baryElem; coordv(7,:); coordv(1,:)]; % [3*singleNE x 2]
    elseif singleNE == 8
        coordTri0Elem = ...
            [baryElem; coordv(1,:); coordv(2,:);
            baryElem; coordv(2,:); coordv(3,:);
            baryElem; coordv(3,:); coordv(4,:);
            baryElem; coordv(4,:); coordv(5,:);
            baryElem; coordv(5,:); coordv(6,:);
            baryElem; coordv(6,:); coordv(7,:);
            baryElem; coordv(7,:); coordv(8,:);
            baryElem; coordv(8,:); coordv(1,:)]; % [3*singleNE x 2]
    elseif singleNE == 9
        coordTri0Elem = ...
            [baryElem; coordv(1,:); coordv(2,:);
            baryElem; coordv(2,:); coordv(3,:);
            baryElem; coordv(3,:); coordv(4,:);
            baryElem; coordv(4,:); coordv(5,:);
            baryElem; coordv(5,:); coordv(6,:);
            baryElem; coordv(6,:); coordv(7,:);
            baryElem; coordv(7,:); coordv(8,:);
            baryElem; coordv(8,:); coordv(9,:);
            baryElem; coordv(9,:); coordv(1,:)]; % [3*singleNE x 2]
    elseif singleNE == 10
        coordTri0Elem = ...
            [baryElem; coordv(1,:); coordv(2,:);
            baryElem; coordv(2,:); coordv(3,:);
            baryElem; coordv(3,:); coordv(4,:);
            baryElem; coordv(4,:); coordv(5,:);
            baryElem; coordv(5,:); coordv(6,:);
            baryElem; coordv(6,:); coordv(7,:);
            baryElem; coordv(7,:); coordv(8,:);
            baryElem; coordv(8,:); coordv(9,:);
            baryElem; coordv(9,:); coordv(10,:);
            baryElem; coordv(10,:); coordv(1,:)]; % [3*singleNE x 2]
    end % if singleNE == 3
    
else
    if singleNE == 3
        coordTri0Elem = coordv; % [3 x 2]
    elseif singleNE == 4
        coordTri0Elem = ...
            [coordv(m(n),:); coordv(m(n+1),:); coordv(m(n+2),:); ...
            coordv(m(n),:); coordv(m(n+2),:); coordv(m(n+3),:)];
    elseif singleNE == 5
        coordTri0Elem = ...
            [coordv(m(n),:); coordv(m(n+1),:); coordv(m(n+2),:); ...
            coordv(m(n),:); coordv(m(n+2),:); coordv(m(n+3),:); ...
            coordv(m(n),:); coordv(m(n+3),:); coordv(m(n+4),:)];
    elseif singleNE == 6
        coordTri0Elem = ...
            [coordv(m(n),:); coordv(m(n+1),:); coordv(m(n+2),:); ...
            coordv(m(n),:); coordv(m(n+2),:); coordv(m(n+3),:); ...
            coordv(m(n),:); coordv(m(n+3),:); coordv(m(n+4),:); ...
            coordv(m(n),:); coordv(m(n+4),:); coordv(m(n+5),:)];
    elseif singleNE == 7
        coordTri0Elem = ...
            [coordv(m(n),:); coordv(m(n+1),:); coordv(m(n+2),:); ...
            coordv(m(n),:); coordv(m(n+2),:); coordv(m(n+3),:); ...
            coordv(m(n),:); coordv(m(n+3),:); coordv(m(n+4),:); ...
            coordv(m(n),:); coordv(m(n+4),:); coordv(m(n+5),:); ...
            coordv(m(n),:); coordv(m(n+5),:); coordv(m(n+6),:)];  
    elseif singleNE == 8
        coordTri0Elem = ...
            [coordv(m(n),:); coordv(m(n+1),:); coordv(m(n+2),:); ...
            coordv(m(n),:); coordv(m(n+2),:); coordv(m(n+3),:); ...
            coordv(m(n),:); coordv(m(n+3),:); coordv(m(n+4),:); ...
            coordv(m(n),:); coordv(m(n+4),:); coordv(m(n+5),:); ...
            coordv(m(n),:); coordv(m(n+5),:); coordv(m(n+6),:); ...
            coordv(m(n),:); coordv(m(n+6),:); coordv(m(n+7),:)];  
    elseif singleNE == 9
        coordTri0Elem = ...
            [coordv(m(n),:); coordv(m(n+1),:); coordv(m(n+2),:); ...
            coordv(m(n),:); coordv(m(n+2),:); coordv(m(n+3),:); ...
            coordv(m(n),:); coordv(m(n+3),:); coordv(m(n+4),:); ...
            coordv(m(n),:); coordv(m(n+4),:); coordv(m(n+5),:); ...
            coordv(m(n),:); coordv(m(n+5),:); coordv(m(n+6),:); ...
            coordv(m(n),:); coordv(m(n+6),:); coordv(m(n+7),:); ...
            coordv(m(n),:); coordv(m(n+7),:); coordv(m(n+8),:)];
    elseif singleNE == 10
        coordTri0Elem = ...
            [coordv(m(n),:); coordv(m(n+1),:); coordv(m(n+2),:); ...
            coordv(m(n),:); coordv(m(n+2),:); coordv(m(n+3),:); ...
            coordv(m(n),:); coordv(m(n+3),:); coordv(m(n+4),:); ...
            coordv(m(n),:); coordv(m(n+4),:); coordv(m(n+5),:); ...
            coordv(m(n),:); coordv(m(n+5),:); coordv(m(n+6),:); ...
            coordv(m(n),:); coordv(m(n+6),:); coordv(m(n+7),:); ...
            coordv(m(n),:); coordv(m(n+7),:); coordv(m(n+8),:); ...
            coordv(m(n),:); coordv(m(n+8),:); coordv(m(n+9),:)];
    end % if singleNE == 3
    
end % if n==0

end % function getcoordTri0Elem
%%<<-- End sub function 1 ---------------------------------------------------------------

%%>> -- Begin sub function 2 -------------------------------------------------------------
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
