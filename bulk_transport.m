function [UD] = bulk_transport(UD_old,dist_coeff,phiD_old,D,I,Aq,Av,B,N,Param,Grid,Scales,dtD)
% authors: Marc Hesse, Jake Jordan
% date: 11 Dec 2019
% Description: takes an explicit timestep for the advection diffusion
% equation for the transport of a bulk quantity (composition, enthalpy,
% tracer) that partitions in to both the fluid and the matrix
%
% Input:
% uD_old
% dist_coeff
% phiD_old
% D
% I
% Aq
% Av
% B
% N
% Param
% Grid
% Scales
% dtD

% Output:
% uD = bulk quantity at new time

phi_c = Scales.phi_c;

phiD_tilde = phiD_old + (1-phi_c*phiD_old).*dist_coeff/phi_c;

% effective two-phase advection matrix
Ae = Aq*spdiags(1./phiD_tilde,0,Grid.N,Grid.N)+phi_c*Av;

Im = I;
Ex = I - dtD*D*Ae;
UD = solve_lbvp(Im,Ex*UD_old,B,Param.g,N);

