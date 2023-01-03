function uD = femTreatDirichletBCPoisson(gD, meshInfo, basesk)
%
%   input:
%       gD, Dirichlet function;
%       meshInfo, the mesh information;
%       basesk, polynomial degree;
%
%
%	YcZhang 6/5/2018
%
%   Last modified 6/5/2018
%
%

%%------- Initial setting -------%%
%--- Dirichlet dofs
DirEdgeIndex = meshInfo.DirichletEdgeIndex;
DirVdofs = meshInfo.edge(DirEdgeIndex,:);
DirVdofs = unique(DirVdofs(:));

%--- get the Dirichlet value uD at DirVdofs
Dir_Vnode = meshInfo.node(DirVdofs,:);
uD_Vdofs = gD(Dir_Vnode(:,1),Dir_Vnode(:,2));

%--- get the Dirichlet value uD at DirFdofs
if basesk >= 2
    oneFdofs = nchoosek((basesk-2)+1,(basesk-2));
    uD_Fdofs = zeros(length(DirEdgeIndex)*oneFdofs,1);
    for CurrEdge = 1:length(DirEdgeIndex)
        edgeIndex = DirEdgeIndex(CurrEdge);
        ePoint1 = meshInfo.node(meshInfo.edge(edgeIndex,1),:); % [1 x 2] the (x-coordinate,y-coordinate) of ePoint1
        ePoint2 = meshInfo.node(meshInfo.edge(edgeIndex,2),:); % [1 x 2] the (x-coordinate,y-coordinate) of ePoint2
        
        x_coord = (1:oneFdofs)*(ePoint1(1) + ePoint2(1))/(oneFdofs+1); % [1 x oneFdofs]
        y_coord = (1:oneFdofs)*(ePoint1(2) + ePoint2(2))/(oneFdofs+1);
        
        uD_Fdofs((CurrEdge-1)*oneFdofs+(1:oneFdofs)) = gD(x_coord',y_coord'); % [oneFdofs x 1]
    end % for
else
   uD_Fdofs = []; 
end

uD = [uD_Vdofs;
    uD_Fdofs];

end % function