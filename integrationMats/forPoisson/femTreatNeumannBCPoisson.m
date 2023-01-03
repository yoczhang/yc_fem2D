function uN = femTreatNeumannBCPoisson(pde, meshInfo, Gaussformula1D, basesk)
%
%   input:
%       gN, Neumann function;
%       meshInfo, the mesh information;
%       Gaussformula1D, the 1d Gauss quadrature formula, [NFpoints x 2],
%           the first column is the 1d coordiantes of all Gauss points,
%           the second column is the weights of all Gauss points;
%       basesk, polynomial degree;
%
%
%	YcZhang 6/5/2018
%
%   Last modified 7/5/2018
%
%
Nelems = meshInfo.Nelems;
Nedges = meshInfo.Nedges;
Nnodes = meshInfo.Nnodes;
NeuEdgeIndex = meshInfo.NeumannEdgeIndex;

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

uN = zeros(totalDofs,1);

for CurrEdge = 1:length(NeuEdgeIndex)
    edgeIndex = NeuEdgeIndex(CurrEdge);
    %--- note that Neumann boundary is also the boundary edge
    elemIndex1 = meshInfo.edge2elem(edgeIndex,1);
    elemIndex2 = meshInfo.edge2elem(edgeIndex,2);
    elemIndex = max(elemIndex1,elemIndex2);
    elem = meshInfo.elem{elemIndex};
    areaElem = meshInfo.areaElem(elemIndex);
    vertices = meshInfo.node(elem,:);
    nu_edge = meshInfo.nuEdges(edgeIndex,:);
    
    mapMat = [vertices(2,1)-vertices(1,1), vertices(3,1)-vertices(1,1);
        vertices(2,2)-vertices(1,2), vertices(3,2)-vertices(1,2)];
     J_det = 2*areaElem;
    
    %--- edge setting
    edge_hE = meshInfo.areaEdge(edgeIndex);
    ePoint1 = meshInfo.node(meshInfo.edge(edgeIndex,1),:); % [1 x 2] the (x-coordinate,y-coordinate) of ePoint1
    ePoint2 = meshInfo.node(meshInfo.edge(edgeIndex,2),:); % [1 x 2] the (x-coordinate,y-coordinate) of ePoint2
        
    %- 1D Gauss Points
    phyGpoints1DX = (ePoint1(1)-ePoint2(1))/2*Gaussformula1D(:,1) + (ePoint1(1)+ePoint2(1))/2;
    phyGpoints1DY = (ePoint1(2)-ePoint2(2))/2*Gaussformula1D(:,1) + (ePoint1(2)+ePoint2(2))/2;
    phyGweights1D = edge_hE*Gaussformula1D(:,2)/2;
    
    refGpoints1DX = 1/J_det*(   mapMat(2,2)*(phyGpoints1DX-vertices(1,1)) ...
        -mapMat(1,2)*(phyGpoints1DY-vertices(1,2))   );
    refGpoints1DY = 1/J_det*(   mapMat(1,1)*(phyGpoints1DY-vertices(1,2)) ...
        -mapMat(2,1)*(phyGpoints1DX-vertices(1,1))   );
    
    [refPb, ~, ~] = femRefBases2D(refGpoints1DX, refGpoints1DY, basesk); % [NFpoints x NTbases]
    
    
    %--- get the integration 
    %- gN value at face-Gauss-points
    value_gN = nu_edge(1)*pde.ux(phyGpoints1DX,phyGpoints1DY) ...
        + nu_edge(2)*pde.uy(phyGpoints1DX,phyGpoints1DY); % [NFpoints x 1]
    
    %- the final value
    [Row, ~] = femGetRowColPoisson(elemIndex, meshInfo, basesk);
    Row_vec = Row(:,1);
    uN_temp = refPb'*(phyGweights1D.*value_gN); % [NTbases x 1]
    uN(Row_vec) = uN(Row_vec) + uN_temp;
    
end % for CurrElem

end % function 