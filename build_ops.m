function [D,G,I]=build_ops(Grid)
% author: Marc Hesse
% date: 09/08/2014, 09/23/2016, 12/31/2017
% description:
% This function computes the discrete divergence and gradient matrices on a
% regular staggered grid using central difference approximations. The
% discrete gradient assumes homogeneous boundary conditions.
% Input:
% Grid = structure containing all pertinent information about the grid.
% Output:
% D = discrete divergence matrix
% G = discrete gradient matrix
% I = identity matrix
%
% Example call:
% >> Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 10;
% >> Grid = build_grid(Grid);
% >> [D,G,I]=build_ops(Grid);

Nx = Grid.Nx; Ny = Grid.Ny; Nz = Grid.Nz; N = Grid.N;

if (Nx>1) && (Ny>1)  % 2D case
    % One dimensinal divergence
    
    
%     Readable implementation
%     % 1D divergence matrices
%     Dx = spdiags([-ones(Nx,1) ones(Nx,1)]/Grid.dx,[0 1],Nx,Nx+1); % 1D div-matrix in x-dir
%     Dy = spdiags([-ones(Ny,1) ones(Ny,1)]/Grid.dy,[0 1],Ny,Ny+1); % 1D div-matrix in y-dir
%     Ix = speye(Nx); Iy = speye(Ny);  % 1D identities in x and y dirs
%     % 2D Tensor-product divergence matrices
%     Dx = kron(Dx,Iy);  % 2D div-matrix in x-dir
%     Dy = kron(Ix,Dy);  % 2D div-matrix in y-dir
%     % Complete 2D divergence
%     D = [Dx Dy];

%     Implementation avoiding intermediate matrices
    D = [kron( spdiags([-ones(Nx,1) ones(Nx,1)]/Grid.dx,[0 1],Nx,Nx+1),speye(Ny) ) ...
        kron( speye(Nx),spdiags([-ones(Ny,1) ones(Ny,1)]/Grid.dy,[0 1],Ny,Ny+1) )];
    
    dof_f_bnd = [Grid.dof_f_xmin; Grid.dof_f_xmax;... % boundary faces
        Grid.dof_f_ymin; Grid.dof_f_ymax];
elseif (Nx>1) && (Ny==1) % 1D x-direction
    D = spdiags([-ones(Nx,1) ones(Nx,1)]/Grid.dx,[0 1],Nx,Nx+1);
    dof_f_bnd = [Grid.dof_f_xmin; Grid.dof_f_xmax];   % boundary faces
elseif (Nx==1) && (Ny>1) % 1D y-direction
    D = spdiags([-ones(Ny,1) ones(Ny,1)]/Grid.dy,[0 1],Ny,Ny+1);
    dof_f_bnd = [Grid.dof_f_ymin; Grid.dof_f_ymax]-Grid.Nfx;   % boundary faces
end



%% Gradient
% Note this is only true in cartesian coordinates!
% For more general coordinate systems it is worth
% assembling G and D seperately.
G = -D';
G(dof_f_bnd,:) = 0;

%% Identity
I = speye(Grid.N);

%% Adjust divergence for different coordinate systems
if strcmp(Grid.geom,'cylindrical_r')
    Rf = spdiags(Grid.xf,0,Nx+1,Nx+1);
    Rcinv = spdiags(1./Grid.xc,0,Nx,Nx);
    D = Rcinv*D*Rf;
elseif strcmp(Grid.geom,'spherical1D') || strcmp(Grid.geom,'spherical_r')
    Rf = spdiags(Grid.xf.^2,0,Grid.Nx+1,Nx+1);
    Rcinv = spdiags(1./(Grid.xc.^2),0,Nx,Nx);
    D = Rcinv*D*Rf;
elseif strcmp(Grid.geom,'cylindrical_rz') % cylindrical coordinates
    % assumes: x-dir is radial direction and y-dir is the cylinder axis
    Rf = spdiags(Grid.xf,0,Nx+1,Nx+1);
    Rcinv = spdiags(1./Grid.xc,0,Nx,Nx);
    D = [kron( Rcinv*spdiags([-ones(Nx,1) ones(Nx,1)]/Grid.dx,[0 1],Nx,Nx+1)*Rf,speye(Ny) ) ...
        kron( speye(Nx),spdiags([-ones(Ny,1) ones(Ny,1)]/Grid.dy,[0 1],Ny,Ny+1) )];
    
    % Implementation below has become obsolete with tensor products
% elseif strcmp(Grid.geom,'cylindrical_rz')
%     % assumes: y-dir is radial direction 
%     %          simplifies the assembly because grid is ordered y-first
%     % The change in geometry goes into 1D matrix before Dy is reassembled
%     Dx = kron( spdiags([-ones(Nx,1) ones(Nx,1)]/Grid.dx,[0 1],Nx,Nx+1),speye(Ny) );
%     Dy1 = spdiags([-ones(Ny,1) ones(Ny,1)]/Grid.dy,[0 1],Ny,Ny+1);
%     
%     Rf = spdiags(Grid.yf,0,Ny+1,Ny+1);
%     Rcinv = spdiags(1./Grid.yc,0,Ny,Ny);
%     Dy = kron( speye(Nx),Rcinv*Dy1*Rf);
% %     Dy = Rcinv*spdiags([-ones(Ny,1) ones(Ny,1)]/Grid.dy,[0 1],Ny,Ny+1)*Rf;
% %     Dy = spblkdiag(Dy,Nx);
%     D = [Dx Dy];
end
