function [Row, Col] = femGetRowColPoisson(CurrElem, meshInfo, basesk)
%
%   input:
%       Nelems, the number of elems;
%       elemIndx, the scalar-type, the elem index, [1 x 1];
%       elem2edge, the vector-type, the all edges index of elem, [1 x singleNE];
%       NTbases, the number of elem T bases;
%       NFbases, the number of edge F bases.
%
%   output:
%       Row, the row of the matrix, [(NTbases + singleNE*NFbases) x (NTbases + singleNE*NFbases)];
%       Col, the column of the matrix, [(NTbases + singleNE*NFbases) x (NTbases + singleNE*NFbases)];
%
%
%	YxQian 7/5/2018
%
%   Last modified 7/5/2018
%
%

elem = meshInfo.elem{CurrElem}; % [1 x singleNE]
elem2edge = meshInfo.elem2edge{CurrElem}; % [1 x singleNE]
NTbases = nchoosek(basesk+2,basesk);
Nnodes = meshInfo.Nnodes;
Nedges = meshInfo.Nedges;

V_row = elem' * ones(1,NTbases);
V_col = ones(NTbases,1) * elem;

% % if basesk >= 2
% %     oneFdofs = nchoosek((basesk-2)+1,(basesk-2));
% %     singleNE = length(elem2edge);
% %     F_row = zeros(singleNE*oneFdofs, NTbases); % [singleNE*oneFdofs x NTbases]
% %     F_col = zeros(NTbases, singleNE*oneFdofs); % [NTbases x singleNE*oneFdofs]
% %     for CurrEdge = 1:singleNE
% %         edgeIndx = elem2edge(CurrEdge);
% %         F_row((CurrEdge-1)*oneFdofs+(1:oneFdofs),:) = ...
% %             Nnodes + (edgeIndx-1)*oneFdofs + (1 : oneFdofs)' * ones(1,NTbases);
% %         F_col(:,(CurrEdge-1)*oneFdofs+(1 : CurrEdge*oneFdofs) ) = ...
% %             Nnodes + (edgeIndx-1)*oneFdofs + ones(NTbases,1) * (1 : oneFdofs);
% %     end
% % else 
% %     F_row = [];
% %     F_col = [];
% % end

if basesk >= 2
    oneFdofs = nchoosek((basesk-2)+1,(basesk-2));
    edgeIndx = elem2edge'; % [singleNE x 1]
    F_row = Nnodes + (edgeIndx-1) * oneFdofs; % (First d.o.f.) of each face.
    F_row = bsxfun(@plus,F_row,1:oneFdofs); % [singleNE x oneFdofs] (d.o.f. for each face) 
    F_row = F_row';
    F_row = F_row(:)*ones(1,NTbases); % [singleNE*oneFdofs x NTbases]
    F_col = F_row'; % [NTbases x singleNE*oneFdofs]
else 
    F_row = [];
    F_col = [];
end

if basesk >= 3
    oneTdofs = nchoosek((basesk-3)+2,(basesk-3));
    T_row = Nnodes + Nedges*oneFdofs + (CurrElem-1)*oneTdofs ...
        + (1:oneTdofs)' * ones(1,NTbases); % [oneTdofs x NTbases]
    T_col = Nnodes + Nedges*oneFdofs + (CurrElem-1)*oneTdofs ...
        + ones(NTbases,1)*(1:oneTdofs); % [NTbases x oneTdofs]
else
    T_row = [];
    T_col = [];
end 


Row = [V_row;
    F_row;
    T_row];
Col = [V_col, F_col, T_col];
end 
%%<<-- End sub function 1 ------------------------------------------------------------