function [numFlux,numFlux_directComp] = hhoNumericalFlux(meshInfo, Gaussformulas, edgeInd, elemInd, ...
    Uh_T, Uh_F, Tbasesk, Fbasesk)
%
%   %----------------------------------------------------------
%       Ref: 2017 (arXiv), An introduction to Hybrid High-Order methods,
%           chapter 3.2.5, Local conservation and flux continuty.
%   %----------------------------------------------------------
%
%   input:
%       meshInfo, the mesh information;
%       edgeInd, the edge indices on which to compute the numerical flux;
%       elemInd, the elem indices corresponding on the edgeInd;
%           i.e., edgeInd(i) lies on the elemInd(i).
%       uhT, the numerical solution on volumn;
%       uhF, the numerical solution on face;
%       Tbasesk, the volumn bases degree;
%       Fbasesk, the face bases degree;
%
%   output:
%
%
%
%   YcZhang 23/1/2018
%
%   Last modified 23/1/2018
%
%

Gaussformula2D = Gaussformulas{1};
Gaussformula1D = Gaussformulas{2};

CoeffOne = @(x,y) 1 + 0.*x;

RTbasesk = Fbasesk + 1;

NTbases = nchoosek(Tbasesk+2,Tbasesk);
NFbases = nchoosek(Fbasesk+1,Fbasesk);
NRTbases = nchoosek(RTbasesk+2,RTbasesk);

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
    
numFlux = 0;
numFlux_directComp = 0;
for nEdge = 1:length(edgeInd)
    globalEdgeInd = edgeInd(nEdge);
    globalElemInd = elemInd(nEdge);
    singleElem = cell2mat(meshInfo.elem(globalElemInd));
    edges_of_elem = meshInfo.elem2edge{globalElemInd};
    
    if isfield(meshInfo,'PermeabilityCoeffs')
        fluxCoeff = meshInfo.PermeabilityCoeffs(globalElemInd);
    else
        fluxCoeff = 1;
    end
    
    elem1 = meshInfo.edge2elem(globalEdgeInd,1);
    local_e1 = meshInfo.edge2elem(globalEdgeInd,3);
    local_e2 = meshInfo.edge2elem(globalEdgeInd,4);
    if elem1 == globalElemInd
        localEdgeInd = local_e1;
    else
        localEdgeInd = local_e2;
    end
    elemNu = meshInfo.nuEdge0Elem{globalElemInd}(:,localEdgeInd);
    
    local_uhT = Uh_T((globalElemInd-1)*NTbases+1:(globalElemInd)*NTbases); % [NTbases x 1]
    local_Fdofs = (edges_of_elem(:)-1)*NFbases;
    local_Fdofs = bsxfun(@plus, local_Fdofs, 1:NFbases);
    local_Fdofs = reshape(local_Fdofs',length(edges_of_elem(:))*NFbases,1); % the all Fdofs on currElem.
    local_uhF = Uh_F(local_Fdofs); % the uh on all faces of currElem.
        
    %% ------- Part 1, get the information about 
    % 1. physical GaussPoints, 
    % 2. different element bases on phy GaussPoints on element.
    %
    singleNE = size(singleElem,2); % get the number of edges on singleElem.
    concavePointElem_ii = meshInfo.concavePointElem(globalElemInd); % [1 x 1], the concavePointElem of ii-the elem.
    coordv = meshInfo.node(singleElem,:); % [singleNE x 2], the vertices (x-coord, y-coord) of element.
    centroidElem = meshInfo.centroidElem(globalElemInd,:); % [1 x 2], the centroid(xing xin) (x-coord, y-coord) of element.
    coordTri0Elem = getcoordTri0Elem(singleNE, concavePointElem_ii, centroidElem, coordv); 
        %> if singleNE==3, coordTri0Elem, [3 x 2]. else coordTri0Elem, [3*singleNE x 2].
    elem_xT = centroidElem(1);
    elem_yT = centroidElem(2); 
    elem_hT = meshInfo.hElem(globalElemInd); 
    
    %% ------- Part 2
    %%------- the mat and vec setting -------%%
    G1 = zeros(NRTbases,NRTbases); %G1 = G1_x + G1_y;
    G3 = zeros(NRTbases,NTbases);
    G4 = zeros(NTbases,NTbases);
    G5 = zeros(NTbases,NRTbases);
    V1 = zeros(NRTbases,1); % the (Pu_T,1)_T vector
    V2 = zeros(NTbases,1); % the (u_T,1)_T vector
    XF = zeros(NRTbases,singleNE*NFbases);
    I1 = eye(NFbases,NFbases); % this 'I' is the Capital of 'i'.
    
    %%------- the integration on little triangles -------%%
    for nt = 1:(size(coordTri0Elem,1)/3)
        coordTri_nt = coordTri0Elem((nt-1)*3+1:nt*3,:); % [3 x 2], (x-coorc,y-coord) of coordTri_nt.
        [phyGpoints2D, phyGweights2D] = getGaussLocalTri(coordTri_nt, Gaussformula2D);

        % get the T-bases values on 2D quad points
        [RT_Pb, RT_Pbx, RT_Pby] = localBases2D(elem_xT, elem_yT, elem_hT, phyGpoints2D(:,1), phyGpoints2D(:,2), RTbasesk);
            %> RT_Pb, RT_Pbx, RT_Pby, [NTpoints x NRTbases].
        [T_Pb, T_Pbx, T_Pby] = localBases2D(elem_xT, elem_yT, elem_hT, phyGpoints2D(:,1), phyGpoints2D(:,2), Tbasesk);
            %> T_Pb, T_Pbx, T_Pby, [NTpoints x NTbases].
        
        %--- the funcValue may be chosen by case.
        valueCoeffOne_2D = CoeffOne(phyGpoints2D(:,1), phyGpoints2D(:,2));
       
        %--- MatElem
        G1 = G1 ...
            + getMatOnElem(valueCoeffOne_2D, RT_Pbx, RT_Pbx, phyGweights2D) ... 
            + getMatOnElem(valueCoeffOne_2D, RT_Pby, RT_Pby, phyGweights2D); % [NRTbases x NRTbases]
        G3 = G3 ...
            + getMatOnElem(valueCoeffOne_2D, T_Pbx, RT_Pbx, phyGweights2D) ... 
            + getMatOnElem(valueCoeffOne_2D, T_Pby, RT_Pby, phyGweights2D); % [NRTbases x NTbases]
        G4 = G4 + getMatOnElem(valueCoeffOne_2D, T_Pb, T_Pb, phyGweights2D); % [NTbases x NTbases]
        G5 = G5 + getMatOnElem(valueCoeffOne_2D, RT_Pb, T_Pb, phyGweights2D); % [NTbases x NRTbases]
        
        V1 = V1 + RT_Pb' * (phyGweights2D.*valueCoeffOne_2D); % [NRTbases x 1]
        V2 = V2 + T_Pb' * (phyGweights2D.*valueCoeffOne_2D); % [NTbases x 1]
    end % for nt
    
    
    %%------- integrations on Faces -------%%
    XT_temp = G3;
    for CurrEdge = 1:singleNE
        %--- edge setting
        edgeIndex = meshInfo.elem2edge{globalElemInd}(CurrEdge);
        ePoint1 = meshInfo.node(meshInfo.edge(edgeIndex,1),:); % [1 x 2] the (x-coordinate,y-coordinate) of ePoint1
        ePoint2 = meshInfo.node(meshInfo.edge(edgeIndex,2),:); % [1 x 2] the (x-coordinate,y-coordinate) of ePoint2
        edge_xE = (ePoint1(1) + ePoint2(1))/2;
        %edge_yE = (ePoint1(2) + ePoint2(2))/2;
        edge_hE = meshInfo.areaEdge(edgeIndex);
        edge_nu = meshInfo.nuEdge0Elem{globalElemInd}(:,CurrEdge); % [1 x 2], the outward unit normal
        
        %--- 1D Gauss Points
        phyGpoints1DX = (ePoint1(1)-ePoint2(1))/2*Gaussformula1D(:,1) + (ePoint1(1)+ePoint2(1))/2; 
        phyGpoints1DY = (ePoint1(2)-ePoint2(2))/2*Gaussformula1D(:,1) + (ePoint1(2)+ePoint2(2))/2; 
        phyGweights1D = edge_hE*Gaussformula1D(:,2)/2;
        
        %--- some values at 1D GaussPoints
        [~, RT_Pbx, RT_Pby] = localBases2D(elem_xT, elem_yT, elem_hT, phyGpoints1DX, phyGpoints1DY, RTbasesk);
            %> RT_Pb, RT_Pbx, RT_Pby, [NFpoints x NRTbases].
        [T_Pb, ~, ~] = localBases2D(elem_xT, elem_yT, elem_hT, phyGpoints1DX, phyGpoints1DY, Tbasesk);
            %> T_Pb, T_Pbx, T_Pby, [NFpoints x NTbases].
        [F_Pb, ~] = localBases1D(edge_xE, edge_hE, Gaussformula1D(:,1), Fbasesk);
            %> F_Pb, F_Pbx, F_Pby, [NFpoints x NFbases].
        valueCoeffOne_1D = CoeffOne(phyGpoints1DX, phyGpoints1DY); % [NFpoints x 1].
        
        %--- integration matrix
        H1_m = edge_nu(1)*getMatOnElem(valueCoeffOne_1D, F_Pb, RT_Pbx, phyGweights1D) ...
            + edge_nu(2)*getMatOnElem(valueCoeffOne_1D, F_Pb, RT_Pby, phyGweights1D); % [NRTbases x NFbases]
        %H2_m = getMatOnElem(valueCoeffOne_1D, T_Pb, F_Pb, phyGweights1D); % [NFbases x NTbases].
        %H3_m = getMatOnElem(valueCoeffOne_1D, RT_Pb, F_Pb, phyGweights1D); % [NFbases x NRTbases].
        %J1_m = getMatOnElem(valueCoeffOne_1D, F_Pb, F_Pb, phyGweights1D); % [NFbases x NFbases].
        K1_m = edge_nu(1)*getMatOnElem(valueCoeffOne_1D, T_Pb, RT_Pbx, phyGweights1D) ...
            + edge_nu(2)*getMatOnElem(valueCoeffOne_1D, T_Pb, RT_Pby, phyGweights1D); % [NRTbases x NTbases]
        
        %--- other matrix
        XT_temp = XT_temp - K1_m; % [NRTbases x NTbases]
        XF(:,(CurrEdge-1)*NFbases+1:CurrEdge*NFbases) = H1_m; % [NRTbases x singleNE*NFbases]
    end
    G2 = G1; G2(1,:) = V1'; % [NRTbases x NRTbases]
    XT = XT_temp; XT(1,:) = V2'; % [NRTbases x NTbases]
    M_T1 = G2\XT; % [NRTbases x NTbases]
    M_T2 = G2\XF;  % [NRTbases x singleNE*NFbases]
    
    %--- get the Potential Reconstructin Operator dofs
    PROdofs = [M_T1,M_T2]*[local_uhT;local_uhF]; % [(NRTbases) x 1]
    
    %% ------- Part 3, the stabilization term
    %%------- to get the stabilizationTerm -------%%
    tau_Fm = 1; % the \tau_{Fm} coefficient of the stabilizationTerm
    J1 = zeros(singleNE*NFbases,singleNE*NFbases);
    stabilizationTerm = zeros(singleNE*NFbases, NTbases+singleNE*NFbases);
    for CurrEdge = 1:singleNE
        I_Fm = zeros(NFbases,singleNE*NFbases);
        %--- edge setting
        edgeIndex = meshInfo.elem2edge{globalElemInd}(CurrEdge);
        ePoint1 = meshInfo.node(meshInfo.edge(edgeIndex,1),:); % [1 x 2] the (x-coordinate,y-coordinate) of ePoint1
        ePoint2 = meshInfo.node(meshInfo.edge(edgeIndex,2),:); % [1 x 2] the (x-coordinate,y-coordinate) of ePoint2
        edge_xE = (ePoint1(1) + ePoint2(1))/2;
        %edge_yE = (ePoint1(2) + ePoint2(2))/2;
        edge_hE = meshInfo.areaEdge(edgeIndex);
        
        %--- 1D Gauss Points
        phyGpoints1DX = (ePoint1(1)-ePoint2(1))/2*Gaussformula1D(:,1) + (ePoint1(1)+ePoint2(1))/2; 
        phyGpoints1DY = (ePoint1(2)-ePoint2(2))/2*Gaussformula1D(:,1) + (ePoint1(2)+ePoint2(2))/2; 
        phyGweights1D = edge_hE*Gaussformula1D(:,2)/2;
        
        %--- some values at 1D GaussPoints
        [RT_Pb, ~, ~] = localBases2D(elem_xT, elem_yT, elem_hT, phyGpoints1DX, phyGpoints1DY, RTbasesk);
            %> RT_Pb, RT_Pbx, RT_Pby, [NFpoints x NRTbases].
        [T_Pb, ~, ~] = localBases2D(elem_xT, elem_yT, elem_hT, phyGpoints1DX, phyGpoints1DY, Tbasesk);
            %> T_Pb, T_Pbx, T_Pby, [NFpoints x NTbases].
        [F_Pb, ~] = localBases1D(edge_xE, edge_hE, Gaussformula1D(:,1), Fbasesk);
            %> F_Pb, F_Pbx, F_Pby, [NFpoints x NFbases].
        valueCoeffOne_1D = CoeffOne(phyGpoints1DX, phyGpoints1DY); % [NFpoints x 1].
        
        %--- integration matrix          
        H2_m = getMatOnElem(valueCoeffOne_1D, T_Pb, F_Pb, phyGweights1D); % [NFbases x NTbases].
        H3_m = getMatOnElem(valueCoeffOne_1D, RT_Pb, F_Pb, phyGweights1D); % [NFbases x NRTbases].
        J1_m = getMatOnElem(valueCoeffOne_1D, F_Pb, F_Pb, phyGweights1D); % [NFbases x NFbases].
        J1((CurrEdge-1)*NFbases+1:(CurrEdge)*NFbases, ...
            (CurrEdge-1)*NFbases+1:(CurrEdge)*NFbases) = J1_m;
        
        %--- other matrix
        I_Fm(:,(CurrEdge-1)*NFbases+1:CurrEdge*NFbases) = I1; % [NFbases x singleNE*NFbases]
        
        %--- get the final integration matrix
        M_Fm_1 = J1_m\(H2_m + (H3_m - H2_m/G4*G5)*M_T1); % [NFbases x NTbases]
        M_Fm_2 = J1_m\(H3_m - H2_m/G4*G5)*M_T2 - I_Fm; % [NFbases x singleNE*NFbases]
        
        stabilizationTerm = stabilizationTerm ...
            + tau_Fm*edge_hE^(-1) * M_Fm_2' * J1_m * [M_Fm_1, M_Fm_2];
            %> [(NTbases + singleNE*NFbases) x (NTbases + singleNE*NFbases)]
    end % for CurrEdge
    
    %--- get the Boundary Residual Operator BRO dofs.
    BROdofs = -J1\stabilizationTerm*[local_uhT;local_uhF]; % [singleNE*NFbases x 1]
    givenEdge_BROdofs = BROdofs((localEdgeInd-1)*NFbases+1:(localEdgeInd)*NFbases);
    
    
    %% ------- Part II, to get the final Numericla Flux on the given edge.
    %%------- on the localEdgeInd edge
    edgeIndex = meshInfo.elem2edge{globalElemInd}(localEdgeInd);
    ePoint1 = meshInfo.node(meshInfo.edge(edgeIndex,1),:); % [1 x 2] the (x-coordinate,y-coordinate) of ePoint1
    ePoint2 = meshInfo.node(meshInfo.edge(edgeIndex,2),:); % [1 x 2] the (x-coordinate,y-coordinate) of ePoint2
    edge_xE = (ePoint1(1) + ePoint2(1))/2;
    %edge_yE = (ePoint1(2) + ePoint2(2))/2;
    edge_hE = meshInfo.areaEdge(edgeIndex);
    
    %--- 1D Gauss Points
    phyGpoints1DX = (ePoint1(1)-ePoint2(1))/2*Gaussformula1D(:,1) + (ePoint1(1)+ePoint2(1))/2;
    phyGpoints1DY = (ePoint1(2)-ePoint2(2))/2*Gaussformula1D(:,1) + (ePoint1(2)+ePoint2(2))/2;
    phyGweights1D = edge_hE*Gaussformula1D(:,2)/2;
    
    [~, RT_Pbx, RT_Pby] = localBases2D(elem_xT, elem_yT, elem_hT, phyGpoints1DX, phyGpoints1DY, RTbasesk);
        %> RT_Pb, RT_Pbx, RT_Pby, [NFpoints x NRTbases].
	[~, T_Pbx, T_Pby] = localBases2D(elem_xT, elem_yT, elem_hT, phyGpoints1DX, phyGpoints1DY, Tbasesk);
        %> RT_Pb, RT_Pbx, RT_Pby, [NFpoints x NRTbases].
    [F_Pb, ~] = localBases1D(edge_xE, edge_hE, Gaussformula1D(:,1), Fbasesk);
        %> F_Pb, F_Pbx, F_Pby, [NFpoints x NFbases].
    %valueCoeffOne_1D = CoeffOne(phyGpoints1DX, phyGpoints1DY); % [NFpoints x 1].
    
    %--- final Numerical Flux
    numFlux = numFlux + ...
        fluxCoeff * ...
        (-(phyGweights1D'*RT_Pbx*elemNu(1)+phyGweights1D'*RT_Pby*elemNu(2)) * (PROdofs)...
        + phyGweights1D'*F_Pb*givenEdge_BROdofs );

%     numFlux_directComp = numFlux_directComp + ...
%         fluxCoeff * ...
%         (-(phyGweights1D'*RT_Pbx*elemNu(1)+phyGweights1D'*RT_Pby*elemNu(2)) * (PROdofs));
    
    numFlux_directComp = numFlux_directComp + ...
        fluxCoeff * ...
        (-(phyGweights1D'*T_Pbx*elemNu(1)+phyGweights1D'*T_Pby*elemNu(2)) * (local_uhT));
    
end % for CurrEdge 




end % function hhoNumericalFlux
