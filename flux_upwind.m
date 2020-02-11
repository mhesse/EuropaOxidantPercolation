function [A] = flux_upwind(q,Grid) % repo (MDOT)
% author: Marc Hesse
% date: 15 April 2015, 8 Nov 2017, 9 Feb 2020
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
    qn = min(q(1:N),0);
    qp = max(q(2:N+1),0);
    A = spdiags([qp,qn],[-1 0],Grid.N+1,Grid.N);
    if strcmp(Grid.periodic,'x-dir') || strcmp(Grid.periodic,'y-dir')
        A(1,Grid.N)   = max(q(1)       ,0);
        A(Grid.N+1,1) = min(q(Grid.N+1),0);
    end
elseif (Nx>1) && (Ny>1) % 2D    
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
    if strcmp(Grid.periodic,'none') 
        % x-matrices
        Axp = kron(spdiags(ones(Nx,1),-1,Nx+1,Nx),speye(Ny));   % 2D x-positive
        Axn = kron(spdiags(ones(Nx,1),0,Nx+1,Nx),speye(Ny));    % 2D x-negative
        
        % y-matrices
        Ayp = kron(speye(Nx),spdiags(ones(Ny,1),-1,Ny+1,Ny));   % 2D y-positive
        Ayn = kron(speye(Nx),spdiags(ones(Ny,1),0,Ny+1,Ny));    % 2D y-negative
    else % Periodic bnd's 
        % 1D matrices
        Axp1 = spdiags(ones(Nx,1),-1,Nx+1,Nx);  % 1D x-poititve
        Axn1 = spdiags(ones(Nx,1),0,Nx+1,Nx);   % 1D x-negative
        Ayp1 = spdiags(ones(Ny,1),-1,Ny+1,Ny);  % 1D y-positive
        Ayn1 = spdiags(ones(Ny,1),0,Ny+1,Ny);   % 1D y-negative
        % Periodic BC's
        if strcmp(Grid.periodic,'x-dir') || strcmp(Grid.periodic,'xy-dir')
            Axp1(1,Grid.Nx)   = 1;
            Axn1(Grid.Nx+1,1) = 1;
        end
        if strcmp(Grid.periodic,'y-dir') || strcmp(Grid.periodic,'xy-dir')
            Ayp1(1,Grid.Ny)   = 1;
            Ayn1(Grid.Ny+1,1) = 1;
        end
        % 2D matrices
        Axp = kron(Axp1,speye(Ny));                 % 2D x-positive
        Axn = kron(Axn1,speye(Ny));                 % 2D x-negative
        Ayp = kron(speye(Nx),Ayp1);                 % 2D y-positive
        Ayn = kron(speye(Nx),Ayn1);                 % 2D y-negative
    end
    
    
    A = spdiags(max(q,0),0,Nf,Nf)*[Axp; Ayp] + spdiags(min(q,0),0,Nf,Nf)*[Axn; Ayn];
    
end