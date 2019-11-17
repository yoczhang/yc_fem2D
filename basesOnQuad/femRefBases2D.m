function [Pb, Pbx, Pby] = femRefBases2D(x, y, k)
%
%   In this function we compute the values of 2D local bases at reference Gauss-Points.
%
%   We let Nbases denote the number of bases, Npoints denote the number of
%   Gauss-Points.
%
%   input:
%       x, the x-coordinates of reference local Gauss-Points, size: 
%                               a row vector, [Npoints x 1].
%       y, the y-coordinates of reference local Gauss-Points, size: 
%                               a row vector, [Npoints x 1].
%       k, the degree of the polynomial, size: a scalar.
%
%   output:
%       Pb, the value of bases at given coordinates, size: 
%                               a matrix, [Npoints x Nbases].
%       Pbx, the derivative to x of the bases, size:
%                               a matrix, [Npoints x Nbases].
%       Pby, the derivative to y of the bases, size:
%                               a matrix, [Npoints x Nbases].
%              
%
%	YxQian 6/5/2018
%
%   Last modified 6/5/2018
%

[r, c] = size(x);
if r < c % we need to guarantee the input x is a column vector.
    x = x'; 
end
[r, c] = size(y);
if r < c % we need to guarantee the input y is a column vector.
    y = y';
end

vone = ones(length(x),1);
if k == 0
    Pb = 1*vone;
    Pbx = 0*vone;
    Pby = 0*vone;
elseif k == 1
    Pb = [1-x-y, x, y];
    Pbx = [-1*vone, 1*vone, 0*vone];
    Pby = [-1*vone, 0*vone, 1*vone];
elseif k == 2
    Pb = [(1-x-y).*(2.*(1-x-y)-1),...
        x .* (2 .* x - 1), ...
        y .* (2 .* y - 1), ...
        4 .* (1 - x - y) .* x, ... % horizontal edge, shui-ping-bian.
        4 .* x .* y, ... % slope edge, xie-bian.
        4 .* (1 - x - y) .* y]; % vertical edge, chui-zhi-bian.
    Pbx = [4.*x + 4.*y - 3, ...
        4.*x - 1*vone, ...
        0*vone, ...
        4*vone - 4.*y - 8.*x, ... % horizontal edge, shui-ping-bian.
        4.*y, ... % slope edge, xie-bian.
        -4.*y]; % vertical edge, chui-zhi-bian.
    Pby = [4.*x + 4.*y - 3*vone, ...
        0*vone, ...
        4.*y - 1*vone, ... 
        -4.*x, ... % horizontal edge, shui-ping-bian.
        4.*x, ... % slope edge, xie-bian.
        4*vone - 8.*y - 4.*x]; % vertical edge, chui-zhi-bian.
end








end % function
