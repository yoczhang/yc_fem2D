function pde = poissonData(coeff_case, u_case)
%
%   poissonData data for poisson equation
%
%
%	YcZhang 6/5/2018
%
%   Last modified 6/5/2018
%

switch coeff_case
    case 0
        K_11 = 1; K_12 = 0;
        K_21 = 0; K_22 = 1;
        K = [K_11, K_12;
            K_21, K_22];
    case 1
        K_11 = 1; K_12 = 0;
        K_21 = 0; K_22 = 1;
        K = [K_11, K_12;
            K_21, K_22];
end

switch u_case
    case 1
        u = @(x,y) 2*x+3*y;
        ux = @(x,y) 2+0.*x;
        uy = @(x,y) 3+0.*x;
        uxy = @(x,y) 0.*x;
        uxx = @(x,y) 0.*x;
        uyy = @(x,y) 0.*x;
    case 2
        u = @(x,y) x.^2+y.*x;
        ux = @(x,y) 2*x + y;
        uy = @(x,y) x + 0.*x;
        uxy = @(x,y) 1 + 0.*x;
        uxx = @(x,y) 2+0.*x;
        uyy = @(x,y) 0.*x;
    case 3
        u = @(x,y) 2*x.^2.*y;
        ux = @(x,y) 4*x.*y;
        uy = @(x,y) 2*x.^2;
        uxy = @(x,y) 4.*x;
        uxx = @(x,y) 4*y;
        uyy = @(x,y) 0.*x;
    case 4
        u = @(x,y) sin(x.*y);
        ux = @(x,y) y.*cos(x.*y);
        uy = @(x,y) x.*cos(x.*y);
        uxy = @(x,y) cos(x.*y) - x.*y.*sin(x.*y);
        uxx = @(x,y) -y.^2.*sin(x.*y);
        uyy = @(x,y) -x.^2.*sin(x.*y);
    case 5 % domain [0,1]x[0,1], zero Dirichlet B.C.
        u = @(x,y) 10*x.^2.*(1-x).^2.*y.^2.*(1-y).^2;
        ux = @(x,y) 20.*x.*y.^2.*(x - 1).^2.*(y - 1).^2 + 10.*x.^2.*y.^2.*(2.*x - 2).*(y - 1).^2;
        uy = @(x,y) 20.*x.^2.*y.*(x - 1).^2.*(y - 1).^2 + 10.*x.^2.*y.^2.*(2.*y - 2).*(x - 1).^2;
        uxx = @(x,y) 20.*x.^2.*y.^2.*(y - 1).^2 + 20.*y.^2.*(x - 1).^2.*(y - 1).^2 ...
            + 40.*x.*y.^2.*(2.*x - 2).*(y - 1).^2;
        uxy = @(x,y) 20.*x.*y.^2.*(2.*y - 2).*(x - 1).^2 + 20.*x.^2.*y.*(2.*x - 2).*(y - 1).^2 ...
            + 10.*x.^2.*y.^2.*(2.*x - 2).*(2.*y - 2) + 40.*x.*y.*(x - 1).^2.*(y - 1).^2;
        uyy = @(x,y) 20.*x.^2.*y.^2.*(x - 1).^2 + 20.*x.^2.*(x - 1).^2.*(y - 1).^2 ...
            + 40.*x.^2.*y.*(2.*y - 2).*(x - 1).^2;
    case 6
        u = @(x,y) exp(-x-y.^2);
        ux = @(x,y) -exp(- y.^2 - x);
        uy = @(x,y) -2*y.*exp(- y.^2 - x);
        uxy = @(x,y) 2.*y.*exp(- y.^2 - x);
        uxx = @(x,y) exp(- y.^2 - x);
        uyy = @(x,y) 4*y.^2.*exp(- y.^2 - x) - 2*exp(- y.^2 - x);
    case 7 % domain [0,1]x[0,1], zero Dirichlet B.C.
        u = @(x,y) (1-y).^2.*(1-x).^2.*sin(x.*y);
        ux = @(x,y) sin(x.*y).*(2.*x - 2).*(y - 1).^2 + y.*cos(x.*y).*(x - 1).^2.*(y - 1).^2;
        uy = @(x,y) sin(x.*y).*(2.*y - 2).*(x - 1).^2 + x.*cos(x.*y).*(x - 1).^2.*(y - 1).^2;
        uxy = @(x,y) cos(x.*y).*(x - 1).^2.*(y - 1).^2 + sin(x.*y).*(2.*x - 2).*(2.*y - 2) + x.*cos(x.*y).*(2.*x - 2).*(y - 1).^2 + y.*cos(x.*y).*(2.*y - 2).*(x - 1).^2 - x.*y.*sin(x.*y).*(x - 1).^2.*(y - 1).^2;
        uxx = @(x,y) 2.*sin(x.*y).*(y - 1).^2 - y.^2.*sin(x.*y).*(x - 1).^2.*(y - 1).^2 + 2.*y.*cos(x.*y).*(2.*x - 2).*(y - 1).^2;
        uyy = @(x,y) 2.*sin(x.*y).*(x - 1).^2 - x.^2.*sin(x.*y).*(x - 1).^2.*(y - 1).^2 + 2.*x.*cos(x.*y).*(2.*y - 2).*(x - 1).^2;
    case 8 % domain [0,1]x[0,1], zero Dirichlet B.C.
        u = @(x,y) sin(pi.*x).^2.*y.*(1+cos(pi.*y)).^2;
        ux = @(x,y) 2.*pi.*y.*cos(pi.*x).*sin(pi.*x).*(cos(pi.*y) + 1).^2;
        uy = @(x,y) sin(pi.*x).^2.*(cos(pi.*y) + 1).^2 - 2.*pi.*y.*sin(pi.*x).^2.*sin(pi.*y).*(cos(pi.*y) + 1);
        uxy = @(x,y) 2.*pi.*cos(pi.*x).*sin(pi.*x).*(cos(pi.*y) + 1).^2 - 4.*pi.^2.*y.*cos(pi.*x).*sin(pi.*x).*sin(pi.*y).*(cos(pi.*y) + 1);
        uxx = @(x,y) 2.*pi.^2.*y.*cos(pi.*x).^2.*(cos(pi.*y) + 1).^2 - 2.*pi.^2.*y.*sin(pi.*x).^2.*(cos(pi.*y) + 1).^2;
        uyy = @(x,y) 2.*pi.^2.*y.*sin(pi.*x).^2.*sin(pi.*y).^2 - 4.*pi.*sin(pi.*x).^2.*sin(pi.*y).*(cos(pi.*y) + 1) - 2.*pi.^2.*y.*cos(pi.*y).*sin(pi.*x).^2.*(cos(pi.*y) + 1);
    case 9
        u = @(x,y) sin(pi.*x).*sin(pi.*y);
        ux = @(x,y) pi.*cos(pi.*x).*sin(pi.*y);
        uy = @(x,y) pi.*cos(pi.*y).*sin(pi.*x);
        uxy = @(x,y) pi.^2.*cos(pi.*x).*cos(pi.*y);
        uxx = @(x,y) -pi.^2.*sin(pi.*x).*sin(pi.*y);
        uyy = @(x,y) -pi.^2.*sin(pi.*x).*sin(pi.*y);
    case 10
        % this case example is from: 
        % LiRui, 2017 A Weak Galerkin Finite Element Method for a Coupled Stokes-Darcy Problem
        % example 1.
        % StokesDomain:[0,1] x [1,2]; DarcyDomain:[0,1] x [0,1].
        u = @(x,y) (2-pi*sin(pi*x)).*(1-y-cos(pi*y));
        ux = @(x,y) pi.^2.*cos(pi.*x).*(y + cos(pi.*y) - 1);
        uy = @(x,y) -(pi.*sin(pi.*x) - 2).*(pi.*sin(pi.*y) - 1);
        uxy = @(x,y) -pi.^2.*cos(pi.*x).*(pi.*sin(pi.*y) - 1);
        uxx = @(x,y) -pi.^3.*sin(pi.*x).*(y + cos(pi.*y) - 1); 
        uyy = @(x,y) -pi.^2.*cos(pi.*y).*(pi.*sin(pi.*x) - 2);
        
    case 11
        % this case example is from: 
        % LiRui, 2017 A Weak Galerkin Finite Element Method for a Coupled Stokes-Darcy Problem
        % example 2.
        % StokesDomain:[0,1] x [1,2]; DarcyDomain:[0,1] x [0,1].
        u = @(x,y) -pi/4*cos(pi*x/2).*y;
        ux = @(x,y) (y.*pi.^2.*sin((pi.*x)/2))/8;
        uy = @(x,y) -(pi.*cos((pi.*x)/2))/4;
        uxy = @(x,y) (pi.^2.*sin((pi.*x)/2))/8;
        uxx = @(x,y) (y.*pi.^3.*cos((pi.*x)/2))/16; 
        uyy = @(x,y) 0.*x;
        
    case 12
        u = @(x,y) x.*y.*(1-x/2).*(1-y/2).*exp(x+y);
        ux = @(x,y) (x.*y.*exp(x + y).*(y/2 - 1))/2 + y.*exp(x + y).*(x/2 - 1).*(y/2 - 1) + x.*y.*exp(x + y).*(x/2 - 1).*(y/2 - 1);
        uy = @(x,y) (x.*y.*exp(x + y).*(x/2 - 1))/2 + x.*exp(x + y).*(x/2 - 1).*(y/2 - 1) + x.*y.*exp(x + y).*(x/2 - 1).*(y/2 - 1);
        uxy = @(x,y) (x.*y.*exp(x + y))/4 + exp(x + y).*(x/2 - 1).*(y/2 - 1) + (x.*exp(x + y).*(y/2 - 1))/2 + (y.*exp(x + y).*(x/2 - 1))/2 + (x.*y.*exp(x + y).*(x/2 - 1))/2 + (x.*y.*exp(x + y).*(y/2 - 1))/2 + x.*exp(x + y).*(x/2 - 1).*(y/2 - 1) + y.*exp(x + y).*(x/2 - 1).*(y/2 - 1) + x.*y.*exp(x + y).*(x/2 - 1).*(y/2 - 1);
        uxx = @(x,y) y.*exp(x + y).*(y/2 - 1) + x.*y.*exp(x + y).*(y/2 - 1) + 2.*y.*exp(x + y).*(x/2 - 1).*(y/2 - 1) + x.*y.*exp(x + y).*(x/2 - 1).*(y/2 - 1);
        uyy = @(x,y) x.*exp(x + y).*(x/2 - 1) + x.*y.*exp(x + y).*(x/2 - 1) + 2.*x.*exp(x + y).*(x/2 - 1).*(y/2 - 1) + x.*y.*exp(x + y).*(x/2 - 1).*(y/2 - 1);
end


f = @(x,y) -(K_11*uxx(x,y)+K_12*uxy(x,y) + K_21*uxy(x,y)+K_22*uyy(x,y));

gD = u;
%gD = @(x,y) 1+ 0.*x;

funcZero = @(x,y) 0.*x;
funcOne = @(x,y) 1+0.*x;

pde = struct(...
    'K_11',K_11, 'K_12', K_12, ...
    'K_21',K_21, 'K_22', K_22, ...
    'K',K, ...
    'u',u, ...
    'ux',ux, ...
    'uy',uy, ...
    'f',f, ...
    'gD',gD, ... % Dirichlet function
    'funcZero',funcZero, ...
    'funcOne',funcOne ...
    );

end % pde function