function [L2_Ferr, H1_Ferr] = hhoL2H1_Ferror(P,Ph,meshInfo,Gaussformula1D,Fbasesk)
%
%
%
%
%
%
%
%
%
%   YcZhang 9/12/2017
%
%   Last modified 9/12/2017
%

NFbases = nchoosek(Fbasesk+1,Fbasesk);

node = meshInfo.node;
edge = meshInfo.edge;
areaEdge = meshInfo.areaEdge;
Nedges = meshInfo.Nedges;

refGpoints = Gaussformula1D(:,1); % [NFpoints x 1]
phyGpoints1DX = (node(edge(:,1),1)-node(edge(:,2),1))/2*refGpoints' ...
    + (node(edge(:,1),1)+node(edge(:,2),1))/2*ones(size(refGpoints')); % [Nedges x NFpoints]
phyGpoints1DY = (node(edge(:,1),2)-node(edge(:,2),2))/2*refGpoints' ...
    + (node(edge(:,1),2)+node(edge(:,2),2))/2*ones(size(refGpoints')); % [Nedges x NFpoints]
phyGpoints1DX = phyGpoints1DX'; % [NFpoints x Nedges]
phyGpoints1DY = phyGpoints1DY';  % [NFpoints x Nedges]
phyGweights1D = areaEdge*(Gaussformula1D(:,2)/2)'; % [Nedges x NFpoints]

edge_xE = 0; edge_hE = 0;
[FPb, FPbx] = localBases1D(edge_xE, edge_hE, refGpoints, Fbasesk); % [NFpoints x NFbases]
wFPb = bsxfun(@times, Gaussformula1D(:,2), FPb); % [NFpoints x NFbases]
wFPbx = bsxfun(@times, Gaussformula1D(:,2), FPbx); % [NFpoints x NFbases]

Pt_h = ((wFPb'*FPb)\wFPb')*P(phyGpoints1DX, phyGpoints1DY); % [NFbases x Nedges]
Ph = reshape(Ph,NFbases,Nedges); % [NFbases x Nedges]

L2_Ferr_temp = phyGweights1D*(FPb*Pt_h - FPb*Ph).^2; % [Nedges x Nedges]
L2_Ferr =  sqrt(sum(diag(L2_Ferr_temp)));

H1_Ferr_temp = phyGweights1D*(FPbx*Pt_h - FPbx*Ph).^2; % [Nedges x Nedges]
H1_Ferr =  sqrt(sum(diag(H1_Ferr_temp)));

end % function