function [D,G,C,I,M]=build_ops(Grid) % repo (MDOT)
% author: Marc Hesse
% date: 09/08/2014, 09/23/2016, 12/31/2017, 01/23/2022
% description:
% This function computes the discrete divergence and gradient matrices on a
% regular staggered grid using central difference approximations. The
% discrete gradient assumes homogeneous boundary conditions.
% Input:
% Grid = structure containing all pertinent information about the grid.
% Output:
% D = discrete divergence matrix
% G = discrete gradient matrix
% C = discrete curl matrix
% I = identity matrix
% M = mean matrix

% Example call:
% >> Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 10;
% >> Grid = build_grid(Grid);
% >> [D,G,C,I,M]=build_ops(Grid);

% Extract from structure for cleaner code
Nx = Grid.Nx; Ny = Grid.Ny; Nz = Grid.Nz; N = Grid.N;
Nfx = Grid.Nfx; Nfy = Grid.Nfy; Nfz = Grid.Nfz; Nf = Grid.Nf;

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

% Natural BC's
if strcmp(Grid.periodic,'none')
    G = -D';
    G(dof_f_bnd,:) = 0;
end

% Periodic BC's
if strcmp(Grid.periodic,'x-dir') && (Nx>1) && (Ny==1)  % periodic in x-direction 1D
    G = spdiags([-ones(Nx,1) ones(Nx,1)]/Grid.dx,[-1 0],Nx+1,Nx);
    G(1,Nx) = -1/Grid.dx; G(Nx+1,1) = 1/Grid.dx;
elseif strcmp(Grid.periodic,'y-dir') && (Ny>1) && (Nx==1)  % periodic in y-direction 1D
    G = spdiags([-ones(Ny,1) ones(Ny,1)]/Grid.dy,[-1 0],Ny+1,Ny);
    G(1,Ny) = -1/Grid.dy; G(Ny+1,1) = 1/Grid.dy;
elseif strcmp(Grid.periodic,'x-dir') && (Nx>1) && (Ny>1) % periodic in x-direction 2D
    Gx = spdiags([-ones(Nx,1) ones(Nx,1)]/Grid.dx,[-1 0],Nx+1,Nx); % 1D grad-matrix in x-dir
    Gx(1,Nx) = -1/Grid.dx; Gx(Nx+1,1) = 1/Grid.dx; % periodic BC's
    Gy = spdiags([-ones(Ny,1) ones(Ny,1)]/Grid.dy,[-1 0],Ny+1,Ny); % 1D grad-matrix in y-dir
    G(1,1) = 0; G(Ny+1,Ny) = 0;    % natural BC's
    Ix = speye(Nx); Iy = speye(Ny);  % 1D identities in x and y dirs
    % 2D Tensor-product divergence matrices
    Gx = kron(Gx,Iy);  % 2D grad-matrix in x-dir
    Gy = kron(Ix,Gy);  % 2D grad-matrix in y-dir
    % Complete 2D divergence
    G = [Gx; Gy];
elseif strcmp(Grid.periodic,'y-dir') && (Nx>1) && (Ny>1) % periodic in y-direction 2D
    Gx = spdiags([-ones(Nx,1) ones(Nx,1)]/Grid.dx,[-1 0],Nx+1,Nx); % 1D grad-matrix in x-dir
    Gx(1,1) = 0; Gx(Nx+1,Nx) = 0; % natural BC's
    Gy = spdiags([-ones(Ny,1) ones(Ny,1)]/Grid.dy,[-1 0],Ny+1,Ny); % 1D grad-matrix in y-dir
    Gy(1,Ny) = -1/Grid.dy; Gy(Ny+1,1) = 1/Grid.dy;   % periodic BC's
    Ix = speye(Nx); Iy = speye(Ny);  % 1D identities in x and y dirs
    % 2D Tensor-product divergence matrices
    Gx = kron(Gx,Iy);  % 2D grad-matrix in x-dir
    Gy = kron(Ix,Gy);  % 2D grad-matrix in y-dir
    % Complete 2D divergence
    G = [Gx; Gy];
elseif strcmp(Grid.periodic,'xy-dir') && (Nx>1) && (Ny>1) % periodic in both directions 2D
    Gx = spdiags([-ones(Nx,1) ones(Nx,1)]/Grid.dx,[-1 0],Nx+1,Nx); % 1D grad-matrix in x-dir
    Gx(1,Nx) = -1/Grid.dx; Gx(Nx+1,1) = 1/Grid.dx; % periodic BC's
    Gy = spdiags([-ones(Ny,1) ones(Ny,1)]/Grid.dy,[-1 0],Ny+1,Ny); % 1D grad-matrix in y-dir
    Gy(1,Ny) = -1/Grid.dy; Gy(Ny+1,1) = 1/Grid.dy;   % periodic BC's
    Ix = speye(Nx); Iy = speye(Ny);  % 1D identities in x and y dirs
    % 2D Tensor-product divergence matrices
    Gx = kron(Gx,Iy);  % 2D grad-matrix in x-dir
    Gy = kron(Ix,Gy);  % 2D grad-matrix in y-dir
    % Complete 2D divergence
    G = [Gx; Gy];
end

%% Identity
I = speye(Grid.N);

%% Curl 
% (currently not implemented)
C = nan;

%% Mean 
if (Nx>1) && (Ny>1)  % 2D case
    Mx = spdiags([ones(Nx,1) ones(Nx,1)]/2,[-1 0],Nx+1,Nx); % 1D mean-matrix in x-dir
    My = spdiags([ones(Ny,1) ones(Ny,1)]/2,[-1 0],Ny+1,Ny); % 1D mean-matrix in y-dir
    Mx = kron(Mx,speye(Ny));                                % 2D mean-matrix in x-dir
    My = kron(speye(Nx),My);                                % 2D mean-matrix in y-dir
    M = [Mx;My];                                            % 2D mean-matrix
    dof_bnd = [Grid.dof_xmin;Grid.dof_xmax;...        % 0 on bnd's
                 Grid.dof_f_ymin;Grid.dof_f_ymax];
    
    dof_f_bnd = [Grid.dof_xmin;Grid.dof_xmax;...        % 0 on bnd's
                 Grid.dof_ymin;Grid.dof_ymax];
   
elseif (Nx>1) && (Ny==1) % 1D x-direction
    M = spdiags([ones(Nx,1) ones(Nx,1)]/2,[-1 0],Nx+1,Nx); % 1D mean-matrix in x-dir
    M(1,1) = 1; M(Nfx,Nx) = 1;
    dof_bnd = [Grid.dof_xmin; Grid.dof_xmax];         % boundary cells
    dof_f_bnd = [Grid.dof_f_xmin; Grid.dof_f_xmax];   % boundary faces
elseif (Nx==1) && (Ny>1) % 1D y-direction
    M = spdiags([ones(Ny,1) ones(Ny,1)]/2,[-1 0],Ny+1,Ny);     % 1D mean-matrix in y-dir
    dof_bnd = [Grid.dof_ymin; Grid.dof_ymax]-Grid.Nfx;         % boundary cells
    dof_f_bnd = [Grid.dof_f_ymin; Grid.dof_f_ymax]-Grid.Nfx;   % boundary faces
end

% for spherical shell coordinates we need periodicity in the y-direction
if strcmp(Grid.geom,'spherical_shell_theta_phi')
    Mx = spdiags([ones(Nx,1) ones(Nx,1)]/2,[-1 0],Nx+1,Nx); % 1D mean-matrix in x-dir
    My = spdiags([ones(Ny,1) ones(Ny,1)]/2,[-1 0],Ny+1,Ny); % 1D mean-matrix in y-dir
    My(1,Ny) = 1/2; My(Ny+1,1) = 1/2;                       % periodicity in y-dir
    Mx = kron(Mx,speye(Ny));                                % 2D mean-matrix in x-dir
    My = kron(speye(Nx),My);                                % 2D mean-matrix in y-dir
    M = [Mx;My];                                            % 2D mean-matrix
    dof_bnd = [Grid.dof_xmin;Grid.dof_xmax];                % boundary cells
    dof_f_bnd = [Grid.dof_f_xmin;Grid.dof_f_xmax];          % boundary faces
end
% Need to think how to do this!!!
% M(dof_f_bnd,dof_bnd) = 1;

%% Adjust operators for different coordinate systems
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
elseif strcmp(Grid.geom,'spherical_shell')
    Rf = spdiags(sin(Grid.xf),0,Nx+1,Nx+1);
    Rcinv = spdiags(1./(Grid.R_shell*sin(Grid.xc)),0,Nx,Nx);
    D = Rcinv*D*Rf;
    G = G/Grid.R_shell;
elseif strcmp(Grid.geom,'spherical_shell_theta_phi')
    % assumes x is polar angle and y is the azimutal angle
    Dx = spdiags([-ones(Nx,1) ones(Nx,1)]/Grid.dx,[0 1],Nx,Nx+1);  % 1D div-matrix in x-dir
    Gx = spdiags([-ones(Nx,1) ones(Nx,1)]/Grid.dx,[-1 0],Nx+1,Nx); % 1D grad-matrix in x-dir
    Dy = spdiags([-ones(Ny,1) ones(Ny,1)]/Grid.dy,[0 1],Ny,Ny+1);  % 1D div-matrix in y-dir
    Gy = spdiags([-ones(Ny,1) ones(Ny,1)]/Grid.dy,[-1 0],Ny+1,Ny); % 1D grad-matrix in y-dir
    Gy(1,Ny) = -1/Grid.dy; Gy(Ny+1,1) = 1/Grid.dy;   % periodic BC's in y-dir
    
    % polar angle
    Sin_f = spdiags(sin(Grid.xf),0,Nx+1,Nx+1);
    R_sin_c_inv = spdiags(1./(Grid.R_shell*sin(Grid.xc)),0,Nx,Nx);
    Dx = R_sin_c_inv*Dx*Sin_f;  % 1D polar D
    Gx = Gx/Grid.R_shell;       % 1D polar G
    Dx = kron(Dx,speye(Ny));    % 2D polar D
    Gx = kron(Gx,speye(Ny));    % 2D polar G
    
    % azimuthal angle
    % Note: The the R matrix is a cell centers for both Dy and Gy
    %       because it is the cell center in the x-dir!
    R_sin_c_inv = spdiags(1./(Grid.R_shell*sin(Grid.xc)),0,Nx,Nx);
    Dy = kron(R_sin_c_inv,Dy);    % 2D azimutal D
    Gy = kron(R_sin_c_inv,Gy);    % 2D azimutal G
    D = [Dx,Dy];
    G = [Gx;Gy];
    
    dof_f_bnd = [Grid.dof_f_xmin; Grid.dof_f_xmax]; % Natural BC's x-dir
    G(dof_f_bnd,:) = 0;
% else
%     fprintf('Currently implemented:\n * cylindrical_r\n * spherical1D\n * cylindrical_rz\n * spherical_shell\n')
%     error('Unknown geometry.' )
end
