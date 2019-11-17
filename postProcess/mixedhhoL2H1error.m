function [potentialL2err, fluxL2err, fluxEnergyErr] = mixedhhoL2H1error(pde, ...
    potentialTdofs, fluxTdofs, fluxFdofs, ...
    meshInfo,Gaussformulas,Tbasesk,Fbasesk)
%
%
%   We let Npoints denote the number of Gauss-Points,
%               Nelems denote the number of the elements of Th,
%               NTbases denote the number of LOCAL bases on each K of Th.
%
%   input:
%       meshInfo, the mesh information.
%       formulaGauss2D, the 2d Gauss quadrature formula, size: a matrix, Npoints x 3,
%               the first column is the x-coordinates of all Gauss-Points,
%               the second column is the y-coordinates of all Gauss-Points,
%               the third is the weights of all Gauss-Points.
%       basesType_trial, polynomial degree.
%
%   output:
%       potentialL2err, \| u-uh \|_{L^2(Th)}, a scalar.
%       fluxL2err, \| K \nabla u - K \nabla uh \|_{L^2(Th)}, a scalar.
%
%
%   YcZhang 3/5/2018
%
%   Last modified 3/5/2018
%

%--- pde setting
tensorK = pde.K;

u = pde.u;
ux = pde.ux;
uy = pde.uy;

% Gauss points setting
CoeffOne = @(x,y) 1 + 0.*x;
Gaussformula2D = Gaussformulas{1};
Gaussformula1D = Gaussformulas{2};

% mesh information: interior edges 
nuEdges = meshInfo.nuEdges; % [Nedges x 2]
    %> the unit normal vector of all edges, 
    %> and here, we set e = T1 \intersect T2, and T1_index > T2_index,
    %> then, the unit normal vector of e, is from T1 to T2 (T1-->T2).
Nelems = meshInfo.Nelems;
RTbasesk = Fbasesk + 1;
NTbases = (Tbasesk+1)*(Tbasesk+2)/2;
NRTbases = nchoosek(RTbasesk+2,RTbasesk);
NFbases = nchoosek(Fbasesk+1,Fbasesk);

potentialL2err = 0;
fluxL2err = 0;
fluxEnergyErr = 0;

% the simple build-in function to get the mat
getMatOnElem = @(funcValue, trialValue, testValue, phyGweights) ...
    testValue' * bsxfun(@times, phyGweights.*funcValue, trialValue);
    %> input:
    %>      funcValue, [Npoints x 1], the value the coeffients function at Gauss points.
    %>      trialValue, [Npoints x NTbases_trial], the trial bases values at Gauss points. 
    %>      testValue, [Npoints x NTbases_test], the test bases values at Gauss points. 
    %>
    %> output:
    %>       [NTbases x NTbases], the mat of integration ...
    %
for CurrElem = 1:Nelems
    %% Part I, get the information about 
    % 1. physical GaussPoints, 
    % 2. different element bases on phy GaussPoints on ii-th element.
    %
    %>>-- Begin Part I -------------------------------- DONOT MODIFY ------------------------------
    K_11 = tensorK(1,1); K_12 = tensorK(1,2);
    K_21 = tensorK(2,1); K_22 = tensorK(2,2);
    
    singleElem = meshInfo.elem{CurrElem,1}; % get the ii-th elem from the cell-type structure: elem.
    singleNE = size(singleElem,2); % get the number of edges on singleElem.
    
    concavePointElem_ii = meshInfo.concavePointElem(CurrElem); % [1 x 1], the concavePointElem of ii-the elem.
    coordv = meshInfo.node(singleElem,:); % [singleNE x 2], the vertices (x-coord, y-coord) of element.
    centroidElem = meshInfo.centroidElem(CurrElem,:); % [1 x 2], the centroid(xing xin) (x-coord, y-coord) of element.
    
    coordTri0Elem = getcoordTri0Elem(singleNE, concavePointElem_ii, centroidElem, coordv); % if singleNE==3, coordTri0Elem, [3 x 2]. else coordTri0Elem, [3*singleNE x 2].
    
    elem_xT = centroidElem(1);  elem_yT = centroidElem(2); 
    elem_hT = meshInfo.hElem(CurrElem); 
    
    % on each element
    pot_Tdofs = potentialTdofs((CurrElem-1)*NTbases+1 : CurrElem*NTbases,1); % [NTbases x 1]
    flu_Tdofs = fluxTdofs((CurrElem-1)*(NTbases-1)+1 : CurrElem*(NTbases-1),1); % [NTbases x 1]
    flu_Fdofs = zeros(singleNE*NFbases,1);
    %<<-- End Part I ---------------------------------------------------------------------------------------
    
    %%------- the mat and vec setting -------%%
    fluxG1 = zeros(NRTbases,NRTbases); 
    fluxG2 = zeros(NRTbases,NTbases);
    divG1 = zeros(NTbases,NTbases);
    divG2 = zeros(NTbases,NTbases);
    V1 = zeros(NRTbases,1); % the (Pu_T,1)_T vector
    divXF = zeros(NTbases,singleNE*NFbases);
    fluxXF = zeros(NRTbases,singleNE*NFbases);
    
    u_TMat = zeros(NTbases,NTbases);
    gradu_TMat = zeros(NTbases,NTbases);
    u_TRhs = zeros(NTbases,1);
    gradu_TRhs = zeros(NTbases,1);
    for nt = 1:(size(coordTri0Elem,1)/3)
        coordTri_nt = coordTri0Elem((nt-1)*3+1:nt*3,:); % [3 x 2], (x-coorc,y-coord) of coordTri_nt.
        [phyGpoints2D, phyGweights2D] = getGaussLocalTri(coordTri_nt, Gaussformula2D);

        %--- get the bases values on quad points
        [RTPb, RTPbx, RTPby] = localBases2D(elem_xT, elem_yT, elem_hT, phyGpoints2D(:,1), phyGpoints2D(:,2), RTbasesk);
            %> RT_Pb, RT_Pbx, RT_Pby, [NTpoints x NRTbases].
        [TPb, TPbx, TPby] = localBases2D(elem_xT, elem_yT, elem_hT, phyGpoints2D(:,1), phyGpoints2D(:,2), Tbasesk);
        wTPb = bsxfun(@times, phyGweights2D, TPb); % Pb, Pbx, Pby, [Npoints x NTbases]      
        
        wTPbx = bsxfun(@times, phyGweights2D, TPbx); % [Npoints x (NTbases)]
        wTPby = bsxfun(@times, phyGweights2D, TPby); % [Npoints x (NTbases)]
        
        %--- the funcValue may be chosen by case.
        valueCoeffOne_2D = CoeffOne(phyGpoints2D(:,1), phyGpoints2D(:,2));
        u_GaussValue2D = u(phyGpoints2D(:,1), phyGpoints2D(:,2)); % [Npoints x 1]
        ux_GaussValue2D = ux(phyGpoints2D(:,1), phyGpoints2D(:,2)); % [Npoints x 1]
        uy_GaussValue2D = uy(phyGpoints2D(:,1), phyGpoints2D(:,2)); % [Npoints x 1]
        
        %--- get integration matrices
        %- exact solution
        u_TMat = u_TMat + (wTPb'*TPb); % [NTbases x NTbases]
        u_TRhs = u_TRhs + wTPb'*u_GaussValue2D; % [NTbases x 1]
        
        gradu_TMat = gradu_TMat + ...
            wTPbx'*(K_11*TPbx + K_12*TPby) + ...
            wTPby'*(K_21*TPbx + K_22*TPby); % [(NTbases) x (NTbases)]
        gradu_TRhs = gradu_TRhs + ...
            wTPbx'*(K_11*ux_GaussValue2D + K_12*uy_GaussValue2D) + ...
            wTPby'*(K_21*ux_GaussValue2D + K_22*uy_GaussValue2D); % [(NTbases) x 1]
        
        %- numerical solution
        %- Consistent flux reconstruction
        fluxG1 = fluxG1 ...
            + K_11.*getMatOnElem(valueCoeffOne_2D, RTPbx, RTPbx, phyGweights2D) ... 
            + K_12.*getMatOnElem(valueCoeffOne_2D, RTPby, RTPbx, phyGweights2D) ... 
            + K_21.*getMatOnElem(valueCoeffOne_2D, RTPbx, RTPby, phyGweights2D) ... 
            + K_22.*getMatOnElem(valueCoeffOne_2D, RTPby, RTPby, phyGweights2D); 
            %> [NRTbases x NRTbases]
        fluxG2 = fluxG2 ...
            - getMatOnElem(valueCoeffOne_2D, TPb, RTPb, phyGweights2D); % [NRTbases x NTbases]
        
        %- Discrete divergence
        divG1 = divG1 ...
            + getMatOnElem(valueCoeffOne_2D, TPb, TPb, phyGweights2D); % [NTbases x NTbases]
        divG2 = divG2 - (...
            + K_11.*getMatOnElem(valueCoeffOne_2D, TPbx, TPbx, phyGweights2D) ... 
            + K_12.*getMatOnElem(valueCoeffOne_2D, TPby, TPbx, phyGweights2D) ... 
            + K_21.*getMatOnElem(valueCoeffOne_2D, TPbx, TPby, phyGweights2D) ... 
            + K_22.*getMatOnElem(valueCoeffOne_2D, TPby, TPby, phyGweights2D) ); % [NRTbases x NTbases]
        
        V1 = V1 + RTPb' * (phyGweights2D.*valueCoeffOne_2D); % [NRTbases x 1]
    end % for nt
    
    %--- get the DOFs of exact solution
    %- u
    u_Tdofs= u_TMat\u_TRhs; % [(NTbases) x 1]
        %> the L2 projection of exact solution on Elem.
    %- grad u
    gradu_Tdofs = gradu_TMat(2:end,2:end)\gradu_TRhs(2:end); % [(NTbases-1) x 1]
    gradu_Fdofs = zeros(singleNE*NFbases,1);
    
    for CurrEdge = 1:singleNE
        %--- edge setting
        edgeIndex = meshInfo.elem2edge{CurrElem}(CurrEdge);
        ePoint1 = meshInfo.node(meshInfo.edge(edgeIndex,1),:); % [1 x 2] the (x-coordinate,y-coordinate) of ePoint1
        ePoint2 = meshInfo.node(meshInfo.edge(edgeIndex,2),:); % [1 x 2] the (x-coordinate,y-coordinate) of ePoint2
        edge_xE = (ePoint1(1) + ePoint2(1))/2;
        nuGlobalEdge = nuEdges(edgeIndex,:); % [1 x 2], unit normal vector of Global Edge.
        edge_hE = meshInfo.areaEdge(edgeIndex);
        edge_nu = meshInfo.nuEdge0Elem{CurrElem}(:,CurrEdge); % [1 x 2], the outward unit normal
        
        %--- 1D Gauss Points
        phyGpoints1DX = (ePoint1(1)-ePoint2(1))/2*Gaussformula1D(:,1) + (ePoint1(1)+ePoint2(1))/2; 
        phyGpoints1DY = (ePoint1(2)-ePoint2(2))/2*Gaussformula1D(:,1) + (ePoint1(2)+ePoint2(2))/2; 
        phyGweights1D = edge_hE*Gaussformula1D(:,2)/2;
        
        %--- some values at 1D GaussPoints
        [RTPb, ~, ~] = localBases2D(elem_xT, elem_yT, elem_hT, phyGpoints1DX, phyGpoints1DY, RTbasesk);
            %> RT_Pb, RT_Pbx, RT_Pby, [NFpoints x NRTbases].
        [TPb, ~, ~] = localBases2D(elem_xT, elem_yT, elem_hT, phyGpoints1DX, phyGpoints1DY, Tbasesk);
            %> T_Pb, T_Pbx, T_Pby, [NFpoints x NTbases].
        [FPb, ~] = localBases1D(edge_xE, edge_hE, Gaussformula1D(:,1), Fbasesk);
            %> F_Pb, F_Pbx, F_Pby, [NFpoints x NFbases].
        valueCoeffOne_1D = CoeffOne(phyGpoints1DX, phyGpoints1DY); % [NFpoints x 1].
        
        %--- integration matrix
        divH1_m = (edge_nu(1)*nuGlobalEdge(1) + edge_nu(2)*nuGlobalEdge(2)) ...
            *getMatOnElem(valueCoeffOne_1D, FPb, TPb, phyGweights1D); % [NTbases x NFbases]
        fluxH1_m = (edge_nu(1)*nuGlobalEdge(1) + edge_nu(2)*nuGlobalEdge(2)) ...
            *getMatOnElem(valueCoeffOne_1D, FPb, RTPb, phyGweights1D); % [NRTbases x NFbases]
        
        %--- other matrix
        divXF(:,(CurrEdge-1)*NFbases+1:CurrEdge*NFbases) = divH1_m; % [NTbases x singleNE*NFbases]
        fluxXF(:,(CurrEdge-1)*NFbases+1:CurrEdge*NFbases) = fluxH1_m; % [NRTbases x singleNE*NFbases]
            
        %- exact solution 
        wFPb = bsxfun(@times, phyGweights1D, FPb); % [NFpoints x NFbases]
        ux_GaussValue1D = ux(phyGpoints1DX, phyGpoints1DY); % [NFpoints x 1]
        uy_GaussValue1D = uy(phyGpoints1DX, phyGpoints1DY); % [NFpoints x 1]
        gradu_Fdofs((CurrEdge-1)*NFbases+1:CurrEdge*NFbases) = ...
            (wFPb'*FPb)\( wFPb'*(nuGlobalEdge(1)*(K_11*ux_GaussValue1D + K_12*uy_GaussValue1D) ...
            + nuGlobalEdge(2)*(K_21*ux_GaussValue1D + K_22*uy_GaussValue1D) ) ); % the exact solution P, projection on the face.
        %- numerical solution
        flu_Fdofs((CurrEdge-1)*NFbases+1:CurrEdge*NFbases) = ...
            fluxFdofs((edgeIndex-1)*NFbases+1:edgeIndex*NFbases);
    end
    
    divMT1 = divG1\divG2(:,2:end); % [NTbases x (NTbases-1)]
    divMT2 = divG1\divXF;  % [NTbases x singleNE*NFbases]
    
    fluxG1 = fluxG1(2:end, 2:end);
    fluxG2 = fluxG2(2:end, :);
    fluxXF = fluxXF(2:end, :);
    fluxMT1 = fluxG1\(fluxG2*divMT1);
    fluxMT2 = fluxG1\(fluxG2*divMT2 + fluxXF);
    
    %--- we need the Dofs on CurrElem, presented by RTbasesk basis function.
    gradu2ElemDofs = [fluxMT1, fluxMT2]*[gradu_Tdofs; gradu_Fdofs]; % [(NRTbases-1) x 1]
    flux2ElemDofs = [fluxMT1, fluxMT2]*[flu_Tdofs; flu_Fdofs]; % [(NRTbases-1) x 1]
    
    for nt = 1:(size(coordTri0Elem,1)/3)
        coordTri_nt = coordTri0Elem((nt-1)*3+1:nt*3,:); % [3 x 2], (x-coorc,y-coord) of coordTri_nt.
        [phyGpoints2D, phyGweights2D] = getGaussLocalTri(coordTri_nt, Gaussformula2D);

        % get the bases values on quad points
        [TPb, ~, ~] = localBases2D(elem_xT, elem_yT, elem_hT, ...
            phyGpoints2D(:,1), phyGpoints2D(:,2), Tbasesk);
            %> trialPb, trialPbx, trialPby, [Npoints x NTbases]   
        [~, RTPbx, RTPby] = localBases2D(elem_xT, elem_yT, elem_hT, ...
            phyGpoints2D(:,1), phyGpoints2D(:,2), RTbasesk);
            %> RTPb, RTPbx, RTPby, [Npoints x NRTbases]
        RTPbx = RTPbx(:,2:end); % [Npoints x (NRTbases-1)]
        RTPby = RTPby(:,2:end); % [Npoints x (NRTbases-1)]
        
        %--- potential L2err
        potentialL2err = potentialL2err + phyGweights2D'*(TPb*u_Tdofs - TPb*pot_Tdofs).^2;
        
        %--- flux L2err, part I
        fluxL2err = fluxL2err ...
            + K_11*phyGweights2D'*(RTPbx*gradu2ElemDofs - RTPbx*flux2ElemDofs).^2 ...
            + K_12*phyGweights2D'*( (RTPby*gradu2ElemDofs - RTPby*flux2ElemDofs).*(RTPbx*gradu2ElemDofs - RTPbx*flux2ElemDofs) ) ...
            + K_21*phyGweights2D'*( (RTPbx*gradu2ElemDofs - RTPbx*flux2ElemDofs).*(RTPby*gradu2ElemDofs - RTPby*flux2ElemDofs) ) ...
            + K_22*phyGweights2D'*(RTPby*gradu2ElemDofs - RTPby*flux2ElemDofs).^2;
    end % for nt
    
    for CurrEdge = 1:singleNE
        %--- edge setting
        edgeIndex = meshInfo.elem2edge{CurrElem}(CurrEdge);
        ePoint1 = meshInfo.node(meshInfo.edge(edgeIndex,1),:); % [1 x 2] the (x-coordinate,y-coordinate) of ePoint1
        ePoint2 = meshInfo.node(meshInfo.edge(edgeIndex,2),:); % [1 x 2] the (x-coordinate,y-coordinate) of ePoint2
        edge_xE = (ePoint1(1) + ePoint2(1))/2;
        nuGlobalEdge = nuEdges(edgeIndex,:); % [1 x 2], unit normal vector of Global Edge.
        edge_hE = meshInfo.areaEdge(edgeIndex);
        %edge_nu = meshInfo.nuEdge0Elem{CurrElem}(:,CurrEdge); % [1 x 2], the outward unit normal
        K_nu_nu = (K_11*nuGlobalEdge(1) + K_12*nuGlobalEdge(2))*nuGlobalEdge(1) ...
            + (K_21*nuGlobalEdge(1) + K_22*nuGlobalEdge(2))*nuGlobalEdge(2);
        
        %--- 1D Gauss Points
        phyGpoints1DX = (ePoint1(1)-ePoint2(1))/2*Gaussformula1D(:,1) + (ePoint1(1)+ePoint2(1))/2; 
        phyGpoints1DY = (ePoint1(2)-ePoint2(2))/2*Gaussformula1D(:,1) + (ePoint1(2)+ePoint2(2))/2; 
        phyGweights1D = edge_hE*Gaussformula1D(:,2)/2;
        
        %--- some values at 1D GaussPoints
        [~, RTPbx, RTPby] = localBases2D(elem_xT, elem_yT, elem_hT, phyGpoints1DX, phyGpoints1DY, RTbasesk);
            %> RT_Pb, RT_Pbx, RT_Pby, [NFpoints x NRTbases].
        [FPb, ~] = localBases1D(edge_xE, edge_hE, Gaussformula1D(:,1), Fbasesk);
            %> F_Pb, F_Pbx, F_Pby, [NFpoints x NFbases].
            
        RTPbx = RTPbx(:,2:end); % [NFpoints x (NRTbases-1)]
        RTPby = RTPby(:,2:end); % [NFpoints x (NRTbases-1)]
        gradu_CurrFdofs = gradu_Fdofs((CurrEdge-1)*NFbases+1:CurrEdge*NFbases); % [NFpoints x 1]
        f_CurrFdofs = flu_Fdofs((CurrEdge-1)*NFbases+1:CurrEdge*NFbases); % [NFpoints x 1]
        
        %--- fluxL2err, part II
        fluxL2err = fluxL2err ...
            + edge_hE*K_nu_nu^(-1)*phyGweights1D'*( ...
            nuGlobalEdge(1)*(K_11*(RTPbx*gradu2ElemDofs-RTPbx*flux2ElemDofs) + K_12*(RTPby*gradu2ElemDofs-RTPby*flux2ElemDofs) ) ...
            + nuGlobalEdge(2)*(K_21*(RTPbx*gradu2ElemDofs-RTPbx*flux2ElemDofs) + K_22*(RTPby*gradu2ElemDofs-RTPby*flux2ElemDofs) ) ...
            - (FPb*f_CurrFdofs-FPb*gradu_CurrFdofs) ...
            ).^2;
        
        %--- flux energy norm
        %- part I
        fluxEnergyErr = fluxEnergyErr + fluxL2err;
        %-part II
        fluxEnergyErr = fluxEnergyErr + ...
            edge_hE*phyGweights1D'*( ...
            FPb*f_CurrFdofs-FPb*gradu_CurrFdofs ).^2;
    end

end % for CurrElem

potentialL2err = sqrt(potentialL2err);
fluxL2err = sqrt(fluxL2err);
fluxEnergyErr = sqrt(fluxEnergyErr);

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
