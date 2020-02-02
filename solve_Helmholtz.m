function [hD,pD,qD] = solve_Helmholtz(D,G,I,phi,n,m,Grid,B,N,fn,Param,Zc) 
% author: Marc Hesse
% date: 29 Mar 2019

% Input: (to be completed)
% D = N by Nf discrete divergence operator
% G = Nf by N discrete gradient operator
% I = N by N identity
% phi = Ny by Nx matrix with porosity field
% n = permeability exponent
% m = compaction exponent
% B = Nc by N constraint matrix
% N = N by N-Nc basis for nullspace of constraint matrix
% fn = r.h.s. vector for Neuman BC
% Param = struture with releavant info about problem parameters
% Zc = N by 1 vector containing the vertical coordinate of cell centers

% Output:
% hD = N by 1 column vector of dimensionless overpressure heads
% pD = N by 1 column vector of dimensionless overpressure
% qD = Nf by 1 column vector of dimensionless relative fluid fluxes

% Porosity Matrices
Phi_n = comp_mean(reshape(phi,Grid.Ny,Grid.Nx).^n,-1,1,Grid);
Phi_m = spdiags(phi.^m,0,Grid.N,Grid.N);

%% Solve mod. Helmholtz equations 
L  = -D*Phi_n*G+Phi_m;              % system matrix
fs = Phi_m*Zc(:); % r.h.s.
flux = @(h) -Phi_n*G*h;
res = @(h,cell) L(cell,:)*h - fs(cell); 

%% Solve boundary value problem
hD = solve_lbvp(L,fs+fn,B,Param.g,N);

% Set proper Dirichlet BC's (hack for changing domain!)
Param2.dof_dir      = [];
Param2.dof_f_dir    = [];
Param2.g            = [];
Param2.dof_neu      = [];
Param2.dof_f_neu    = [];
Param2.qb            = [];

qD = comp_flux_gen(flux,res,hD,Grid,Param2);
pD = hD-Zc(:);
end
