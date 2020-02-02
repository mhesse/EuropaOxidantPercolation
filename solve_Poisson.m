function [uD,vD] = solve_Poisson(D,G,I,phi,m,pD,Grid,B,N,fn,Param)
% author: Marc Hesse
% date: 29 Mar 2019

% Input: (to be completed)
% D = N by Nf discrete divergence operator
% G = Nf by N discrete gradient operator
% I = N by N identity
% phi = Ny by Nx matrix with porosity field
% m = compaction exponent
% pD = N by 1 column vector of dimensionless overpressures
% Grid = structure containing useful info about the grid
% B = Nc by N constraint matrix
% N = N by N-Nc basis for nullspace of constraint matrix
% fn = r.h.s. vector for Neuman BC
% Param = struture with releavant info about problem parameters

% Output:
% uD = N  by 1 column vector of dimensionless solid velocity potentials
% vD = Nf by 1 column vector of dimensionless solid velocities

Phi_m = spdiags(phi.^m,0,Grid.N,Grid.N);
L  = -D*G; % system matrix
fs = Phi_m*pD;   % r.h.s.
flux = @(u) -G*u;
res = @(u,cell) L(cell,:)*u - fs(cell);

% Solve boundary value problem
Param2.dof_dir   = [1];
Param2.dof_f_dir = [1];
Param2.g         = [0];
Param2.dof_neu   = [];
Param2.dof_f_neu = [];
Param2.qb        = [];


uD = solve_lbvp(L,fs+fn,B,Param.g,N);
vD = comp_flux_gen(flux,res,uD,Grid,Param2);

end