function get_16mesh_perturbingQuad

close all
clearvars
clc

%%
%------ mesh16
%-------------------- symmetryPolygon Mesh ---------------------
for nn = 1:4
    close all;
    meshN = num2str(2^(nn+1));
    %--- Stokes Domain [0,1]x[0,1]
    meshname = 'perturbingQuad';
    Smeshname = ['Smesh_',meshname,'_[0,1]x[1,2]_',num2str(meshN)];
    load(Smeshname); Snode = node; Snode(:,2) = Snode(:,2)-1; Selem = elem;
    SmeshName = ['Smesh_16mesh_perturbingQuad_[0,1]x[0,1]_',meshN];
    SmeshInfo = polyMeshAuxStructure(Snode, Selem);
    node = SmeshInfo.node; elem = SmeshInfo.elem;
    save(SmeshName, 'node', 'elem')
 
    %--- Darcy Domain [0,1]x[-1,0]
    meshname = 'perturbingQuad'; 
    Dmeshname = ['Dmesh_',meshname,'_[0,1]x[0,1]_',num2str(meshN)];
    load(Dmeshname); Dnode = node; Dnode(:,2) = Dnode(:,2)-1; Delem = elem; 
    DmeshName = ['Dmesh_16mesh_perturbingQuad_[0,1]x[-1,0]_',meshN];
    DmeshInfo = polyMeshAuxStructure(Dnode, Delem);
    node = DmeshInfo.node; elem = DmeshInfo.elem;
    save(DmeshName, 'node', 'elem')
   
    PolyMshr_PlotMsh_twoDomain1(SmeshInfo.node,SmeshInfo.elem,length(Selem),...
        DmeshInfo.node,DmeshInfo.elem,length(Delem))
    figure
    PolyMshr_PlotMsh_twoDomain2(SmeshInfo.node,SmeshInfo.elem,length(Selem),...
        DmeshInfo.node,DmeshInfo.elem,length(Delem))
end % nn



end % function


%%------- sub function
%
%--- sub function 1
function PolyMshr_PlotMsh_twoDomain1(Snode,Selem,SNelem,Dnode,Delem,DNelem)
Sxx = Snode(:,1);
Syy = Snode(:,2);

hold on;
Selem = Selem(1:SNelem)';                 %Only plot the first block
MaxNVer = max(cellfun(@numel,Selem));      %Max. num. of vertices in mesh
PadWNaN = @(E) [E NaN(1,MaxNVer-numel(E))];  %Pad cells with NaN
ElemMat = cellfun(PadWNaN,Selem,'UniformOutput',false);
ElemMat = vertcat(ElemMat{:});               %Create padded element matrix
patch('Faces',ElemMat,'Vertices',Snode,'FaceColor','w'); pause(1e-6)

hold on
Dxx = Dnode(:,1);
Dyy = Dnode(:,2);
Delem = Delem(1:DNelem)';                 %Only plot the first block
MaxNVer = max(cellfun(@numel,Delem));      %Max. num. of vertices in mesh
PadWNaN = @(E) [E NaN(1,MaxNVer-numel(E))];  %Pad cells with NaN
ElemMat = cellfun(PadWNaN,Delem,'UniformOutput',false);
ElemMat = vertcat(ElemMat{:});               %Create padded element matrix
patch('Faces',ElemMat,'Vertices',Dnode,'FaceColor','w'); pause(1e-6)

axis equal; axis([min([Sxx;Dxx])-0.1 max([Sxx;Dxx])+0.1 min([Syy;Dyy])-0.1 max([Syy;Dyy])+0.1])
end


%--- sub function 2
function PolyMshr_PlotMsh_twoDomain2(Snode,Selem,SNelem,Dnode,Delem,DNelem)
Sxx = Snode(:,1);
Syy = Snode(:,2);

hold on;
Selem = Selem(1:SNelem)';                 %Only plot the first block
MaxNVer = max(cellfun(@numel,Selem));      %Max. num. of vertices in mesh
PadWNaN = @(E) [E NaN(1,MaxNVer-numel(E))];  %Pad cells with NaN
ElemMat = cellfun(PadWNaN,Selem,'UniformOutput',false);
ElemMat = vertcat(ElemMat{:});               %Create padded element matrix
patch('Faces',ElemMat,'Vertices',Snode,'FaceColor','w'); pause(1e-6)

hold on
Dxx = Dnode(:,1);
Dyy = Dnode(:,2);
Delem = Delem(1:DNelem)';                 %Only plot the first block
MaxNVer = max(cellfun(@numel,Delem));      %Max. num. of vertices in mesh
PadWNaN = @(E) [E NaN(1,MaxNVer-numel(E))];  %Pad cells with NaN
ElemMat = cellfun(PadWNaN,Delem,'UniformOutput',false);
ElemMat = vertcat(ElemMat{:});               %Create padded element matrix
patch('Faces',ElemMat,'Vertices',Dnode,'FaceColor','w'); pause(1e-6)

axis equal; 
axis off;
%axis([min([Sxx;Dxx])-0.1 max([Sxx;Dxx])+0.1 min([Syy;Dyy])-0.1 max([Syy;Dyy])+0.1])
end


%--- sub function 3
function node = modifiedNode(node,left_bound,right_bound,bottom_bound,top_bound)

xx = node(:,1);
yy = node(:,2);

tol = 6e-9;

left_abs_xx = abs( xx - left_bound );
right_abs_xx = abs( xx - right_bound );
bottom_abs_yy = abs( yy - bottom_bound );
top_abs_yy = abs( yy - top_bound );

node(left_abs_xx<tol,1) = left_bound;
node(right_abs_xx<tol,1) = right_bound;
node(bottom_abs_yy<tol,2) = bottom_bound;
node(top_abs_yy<tol,2) = top_bound;

end % function modifiedNode



