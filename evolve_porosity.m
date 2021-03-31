function [phiD,Av] = evolve_porosity(D,I,phiD,vD,pD,B,N,Param,Grid,Scales,theta,dtD)
% author: Marc Hesse
% date: 16 Apr 2019

% Input:
% D = N by Nf discrete divergence operator
% G = Nf by N discrete gradient operator
% I = N by N identity
% phiD = N by 1 column vector of scaled porosities
% vD = Nf by 1 column vector of dimensionless solid velocities
% pD = N by 1 column vector of dimensionless pressures
% B = Nc by N constraint matrix
% N = N by N-Nc basis for nullspace of constraint matrix
% fn = r.h.s. vector for Neuman BC
% Grid = structure containing info about grid
% Param = struture with releavant info about problem parameters
% Scales = stucture containg all characteristic scalse
% theta = variable determining time integration (1 = FE, .5 = CN, 0 = BE)
% dtD = dimensionless timestep

% Output:
% phiD = N by 1 column vector of dimensionless porosities
% Av = Nf by N upwinf flux matrix

Av = flux_upwind(vD,Grid);
L = Scales.phi_c*D*Av - spdiags(pD,0,Grid.N,Grid.N);
IM = I + dtD*(1-theta)*L;
EX = I - dtD*theta*L;
phiD = solve_lbvp(IM,EX*phiD,B,Param.g,N);
% if max(phiD*Scales.phi_c) > 1 || min(phiD*Scales.phi_c) < 0; error('Porosity outside  [0 1].\n'); end