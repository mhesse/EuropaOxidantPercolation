function [A] = flux_upwind(q,Grid) % repo
% author: Marc Hesse
% date: 15 April 2015, 8 Nov 2017
% Description:
% This function computes the upwind flux matrix from the flux vector.
%
% Input:
% q = Nf by 1 flux vector from the flow problem.
% Grid = structure containing all pertinent information about the grid.
%
% Output:
% A = Nf by Nf matrix contining the upwinded fluxes
%
% Example call:
% >> Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 10;
% >> Grid = build_grid(Grid);
% >> q = ones(Grid.Nf,1);
% >> [A] = flux_upwind(q,Grid);

Nx = Grid.Nx; Ny = Grid.Ny; Nz = Grid.Nz; N = Grid.N;
Nfx = Grid.Nfx; % # of x faces
Nfy = Grid.Nfy; % # of y faces
Nf  = Grid.Nf;   % # faces

if ((Nx>1) && (Ny==1)) || ((Nx==1) && (Ny>1)) % 1D
    %% One dimensinal
    % too make this work for 1D in y-dir need to replace Nx with N!
    qn = min(q(1:N),0);
    qp = max(q(2:N+1),0);
    A = spdiags([qp,qn],[-1 0],Grid.N+1,Grid.N);
elseif (Nx>1) && (Ny>1) % 2D
%     %% Old direct implementation
%     %% Two dimensional
%     % x - fluxes
%     qxn = min(q(1:Nfx-Ny),0); 
%     qxp = max(q(Ny+1:Nfx),0);
%     Ax = spdiags([qxp,qxn],[-Ny 0],Nfx,N);
%     
%     % y-fluxes
%     QY = reshape(q(Nfx+1:end),Grid.Ny+1,Grid.Nx);
%     qyn = min(reshape(QY(1:Grid.Ny,:),Grid.N,1),0);
%     qyp = max(reshape(QY(2:Grid.Ny+1,:),Grid.N,1),0);
%     row_p = (Grid.Ny+1)*repmat([0:Grid.Nx-1],Grid.Ny,1)+repmat([2:Grid.Ny+1]',1,Grid.Nx); 
%     row_n = (Grid.Ny+1)*repmat([0:Grid.Nx-1],Grid.Ny,1)+repmat([1:Grid.Ny]',  1,Grid.Nx); 
%     Ay = sparse([row_p(:);row_n(:)],[Grid.dof;Grid.dof],[qyp;qyn]);
%     
%     A = [Ax; Ay];
    
%     %% Readable implementation with kron
%     % x-matrices
%     Iy = speye(Ny);
%     Axp1 = spdiags(ones(Nx,1),-1,Nx+1,Nx);  % 1D x-poititve
%     Axn1 = spdiags(ones(Nx,1),0,Nx+1,Nx);   % 1D x-negative
%     Axp = kron(Axp1,Iy);                    % 2D x-positive
%     Axn = kron(Axn1,Iy);                    % 2D x-negative
%     
%     % y-matrices
%     Ix = speye(Nx);
%     Ayp1 = spdiags(ones(Ny,1),-1,Ny+1,Ny);  % 1D y-positive
%     Ayn1 = spdiags(ones(Ny,1),0,Ny+1,Ny);   % 1D y-negative
%     Ayp = kron(Ix,Ayp1);                    % 2D y-positive
%     Ayn = kron(Ix,Ayn1);                    % 2D y-negative
%     
%     % Positive and Negative Matrices
%     Ap = [Axp; Ayp];
%     An = [Axn; Ayn];
%     
%     % Diagonal Flux Matrices
%     Qp = spdiags(max(q,0),0,Nf,Nf);
%     Qn = spdiags(min(q,0),0,Nf,Nf);
%     
%     A = Qp*Ap + Qn*An;
   

    %% Compact implementation with kron
    % x-matrices
    Axp = kron(spdiags(ones(Nx,1),-1,Nx+1,Nx),speye(Ny));   % 2D x-positive
    Axn = kron(spdiags(ones(Nx,1),0,Nx+1,Nx),speye(Ny));    % 2D x-negative
    
    % y-matrices
    Ayp = kron(speye(Nx),spdiags(ones(Ny,1),-1,Ny+1,Ny));   % 2D y-positive
    Ayn = kron(speye(Nx),spdiags(ones(Ny,1),0,Ny+1,Ny));    % 2D y-negative
    
    
    A = spdiags(max(q,0),0,Nf,Nf)*[Axp; Ayp] + spdiags(min(q,0),0,Nf,Nf)*[Axn; Ayn];
    
end